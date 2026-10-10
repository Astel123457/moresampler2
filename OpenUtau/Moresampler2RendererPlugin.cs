using System.Runtime.InteropServices;
using NAudio.Wave;
using OpenUtau.Classic;
using OpenUtau.Core;
using OpenUtau.Core.Format;
using OpenUtau.Core.Render;
using OpenUtau.Core.Ustx;
using OpenUtau.Core.Util;
using Serilog;

namespace Moresampler2.OpenUtau;

[ExternalRenderer("io.github.plus-minus-null.moresampler2", "MORESAMPLER2")]
public sealed class Moresampler2RendererPlugin : IOpenUtauRendererPlugin {
    static UExpressionDescriptor Curve(string name, string abbr) =>
        new(name, abbr, -100, 100, 0) { type = UExpressionType.Curve };
    static UExpressionDescriptor Flag(string name, string abbr, string flag, float min = -100, float max = 100) =>
        new(name, abbr, min, max, 0, flag);
    public int ApiVersion => ExternalRendererRegistry.ApiVersion;
    public RendererPluginMetadata Metadata { get; } = new() {
        Capabilities = new RendererCapabilitiesManifest { cancellation = false, parallelism = 1 },
        AnalysisFormats = new Dictionary<string, AnalysisFormatManifest> {
            ["llsm2"] = new() { name = "moresampler2 LLSM Layer-1 analysis", path = "{wav_dir}/{wav_stem}.llsm2", required = true, canGenerate = true },
        },
        Expressions = new Dictionary<string, UExpressionDescriptor> {
            [Ustx.GENC] = Curve("Gender", Ustx.GENC),
            [Ustx.TENC] = Curve("Tension", Ustx.TENC),
            [Ustx.BREC] = Curve("Breathiness", Ustx.BREC),
            [Ustx.VOIC] = new("Voicing", Ustx.VOIC, 0, 100, 100) { type = UExpressionType.Curve },
            ["MSY"] = new("Noise suppression (Y)", "MSY", 0, 100, 100, "Y"),
            ["MDRY"] = Curve("Dryness", "MDRY"),
            ["MOPN"] = Curve("Mouth opening", "MOPN"),
            ["MFRM"] = Curve("Formant emphasis", "MFRM"),
            ["MSFM"] = Curve("Singer formant", "MSFM"),
            ["MGRW"] = Curve("Growl", "MGRW"),
            ["MCRS"] = Curve("Growl coarseness", "MCRS"),
            ["MDST"] = Curve("Distortion", "MDST"),
            ["MPRN"] = Curve("Pitch randomness", "MPRN"),
            ["MUVC"] = Curve("Unvoiced emphasis", "MUVC"),
            ["MAMP"] = Curve("Pitch-linked amplitude", "MAMP"),
            ["MSTA"] = new UExpressionDescriptor("Stabilization", "MSTA", 0, 10, 0) {
                type = UExpressionType.Curve,
            },
        },
        Settings = new Dictionary<string, RendererSettingDescriptor> {
            ["analysisF0Min"] = new() {
                Name = "Minimum detected F0 (Hz)", Description = "Lowest source pitch to detect for new analysis. Existing analysis is reused until manually overwritten. Must be below maximum F0.",
                Type = RendererSettingType.Number, DefaultValue = "50", Min = 20, Max = 1999, Step = 1,
            },
            ["analysisF0Max"] = new() {
                Name = "Maximum detected F0 (Hz)", Description = "Highest source pitch to detect. Must be above minimum F0.",
                Type = RendererSettingType.Number, DefaultValue = "800", Min = 21, Max = 2000, Step = 1,
            },
            ["analysisVoicingBias"] = new() {
                Name = "Analysis voicing bias", Description = "pYIN bias: values above 1 favor voiced frames; values below 1 favor unvoiced frames.",
                Type = RendererSettingType.Number, DefaultValue = "2", Min = .1, Max = 10, Step = .1,
            },
            ["analysisTransitionRange"] = new() {
                Name = "Pitch transition range (bins)", Description = "pYIN transition distance in quantized pitch bins, not semitones. Smaller values discourage large pitch jumps.",
                Type = RendererSettingType.Integer, DefaultValue = "24", Min = 1, Max = 48, Step = 1,
            },
            ["analysisRefineF0"] = new() {
                Name = "Refine detected F0", Description = "Use libllsm2's refinement of the detected source pitch.",
                Type = RendererSettingType.Boolean, DefaultValue = "true",
            },
            ["outputGainDb"] = new() {
                Name = "Output gain (dB)",
                Description = "Adjusts this renderer's phrase output before track effects.",
                Type = RendererSettingType.Number,
                DefaultValue = "0",
                Min = -24,
                Max = 12,
                Step = 0.5,
            },
        },
    };
    public IRenderer CreateRenderer(RendererPluginContext context) => new Moresampler2Renderer(context);
    public IRendererAnalysisProvider CreateAnalysisProvider(RendererPluginContext context) =>
        new Moresampler2AnalysisProvider();
}

sealed class Moresampler2AnalysisProvider : IRendererAnalysisProvider, IDisposable {
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate int AnalyzeNative(IntPtr input, IntPtr output, ref AnalysisOptions options, int overwrite);
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate int ValidateNative(IntPtr analysis);
    readonly IntPtr library;
    readonly AnalyzeNative analyze;
    readonly ValidateNative validate;

    public Moresampler2AnalysisProvider() {
        var name = OperatingSystem.IsWindows() ? "moresampler2.dll" :
            OperatingSystem.IsMacOS() ? "libmoresampler2.dylib" : "libmoresampler2.so";
        library = NativeLibrary.Load(Moresampler2Engine.ExtractNativeLibrary(name));
        RendererProfile.Configure(library);
        analyze = Marshal.GetDelegateForFunctionPointer<AnalyzeNative>(
            NativeLibrary.GetExport(library, "ms2_generate_v1"));
        validate = Marshal.GetDelegateForFunctionPointer<ValidateNative>(
            NativeLibrary.GetExport(library, "ms2_validate_analysis_v1"));
    }

    public Task<IReadOnlyList<RendererAnalysisResult>> GenerateAsync(
            IReadOnlyList<RendererAnalysisRequest> requests,
            IProgress<int> progress, CancellationToken cancellation) =>
        Task.Run<IReadOnlyList<RendererAnalysisResult>>(() => {
            int completed = 0;
            var results = new List<RendererAnalysisResult>(requests.Count);
            foreach (var request in requests) {
                using var timing = RendererProfile.Measure("analysis.generate", request.SourceFile);
                cancellation.ThrowIfCancellationRequested();
                var input = Marshal.StringToCoTaskMemUTF8(request.SourceFile);
                var output = Marshal.StringToCoTaskMemUTF8(request.OutputFile);
                try {
                    var options = AnalysisOptions.From(request.Settings);
                    bool existed = File.Exists(request.OutputFile);
                    bool invalid = existed && !request.Overwrite && validate(output) != 0;
                    if (invalid) Log.Warning("Unusable moresampler2 analysis at {Path}; regenerating from {Source}",
                        request.OutputFile, request.SourceFile);
                    int code = analyze(input, output, ref options, request.Overwrite || invalid ? 1 : 0);
                    RendererProfile.Write("analysis.generate.result", $"code={code} path={request.OutputFile}");
                    results.Add(code == 0
                        ? new RendererAnalysisResult(request, existed && !request.Overwrite && !invalid
                            ? RendererAnalysisOutcome.AlreadyValid : RendererAnalysisOutcome.Generated)
                        : new RendererAnalysisResult(request, RendererAnalysisOutcome.Failed,
                            $"moresampler2 analysis returned {code} for {request.SourceFile}."));
                } catch (Exception exception) {
                    results.Add(new RendererAnalysisResult(
                        request, RendererAnalysisOutcome.Failed, exception.Message));
                } finally {
                    Marshal.FreeCoTaskMem(input);
                    Marshal.FreeCoTaskMem(output);
                }
                progress.Report(++completed);
            }
            return results;
        }, cancellation);

    public ValueTask<RendererAnalysisState> ValidateAsync(
            RendererAnalysisRequest request, CancellationToken cancellation) {
        using var timing = RendererProfile.Measure("analysis.validate", request.OutputFile);
        cancellation.ThrowIfCancellationRequested();
        // Timestamp and settings changes never authorize replacement of existing analysis.
        var basic = File.Exists(request.OutputFile) ? RendererAnalysisState.Valid : RendererAnalysisState.Missing;
        RendererProfile.Write("analysis.basic_state", $"state={basic} path={request.OutputFile}");
        if (basic != RendererAnalysisState.Valid) return ValueTask.FromResult(basic);
        var path = Marshal.StringToCoTaskMemUTF8(request.OutputFile);
        try {
            return ValueTask.FromResult(validate(path) == 0
                ? RendererAnalysisState.Valid : RendererAnalysisState.Invalid);
        } finally {
            Marshal.FreeCoTaskMem(path);
        }
    }

    public void Dispose() => NativeLibrary.Free(library);
}

public sealed class Moresampler2Renderer : IRenderer {
    readonly Moresampler2Engine engine;
    public Moresampler2Renderer(RendererPluginContext context) {
        engine = new Moresampler2Engine(context);
    }
    public USingerType SingerType => USingerType.Classic;
    public bool SupportsRenderPitch => false;
    public bool SupportsRealCurve => false;
    public bool SupportsExpression(UExpressionDescriptor d) => d.isFlag || !string.IsNullOrEmpty(d.flag) ||
        d.abbr is Ustx.DYN or Ustx.PITD or Ustx.VEL or Ustx.VOL or Ustx.MOD or Ustx.GENC or Ustx.TENC or Ustx.BREC or Ustx.VOIC ||
        d.abbr is "MDRY" or "MOPN" or "MFRM" or "MSFM" or "MGRW" or "MCRS" or "MDST" or "MPRN" or "MUVC" or "MAMP" or "MSTA";
    public RenderResult Layout(RenderPhrase phrase) => new() {
        leadingMs = phrase.leadingMs, positionMs = phrase.positionMs,
        estimatedLengthMs = phrase.durationMs + phrase.leadingMs,
    };
    public Task<RenderResult> Render(RenderPhrase phrase, Progress progress, int trackNo,
            CancellationTokenSource cancellation, bool isPreRender = false,
            RenderPhraseEvents? renderEvents = null) => Task.Run(() => {
        cancellation.Token.ThrowIfCancellationRequested();
        var items = phrase.phones.Select(phone => new ResamplerItem(phrase, phone, engine)).ToArray();
        var info = $"Track {trackNo + 1}: moresampler2 \"{string.Join(" ", phrase.phones.Select(p => p.phoneme))}\"";
        progress.Complete(0, info);
        var result = Layout(phrase);
        var synthesisItems = items.Where(item => !item.phone.direct).ToArray();
        result.samples = synthesisItems.Length == 0
            ? new float[(int)Math.Ceiling(result.estimatedLengthMs * 44.1)]
            : engine.RenderPhrase(phrase, synthesisItems, cancellation.Token);
        AddDirects(phrase, items, result.samples);
        if (phrase.rendererSettings.TryGetValue("outputGainDb", out var gainText) &&
                double.TryParse(gainText, System.Globalization.NumberStyles.Float,
                    System.Globalization.CultureInfo.InvariantCulture, out var gainDb)) {
            float gain = (float)Math.Pow(10, Math.Clamp(gainDb, -24, 12) / 20);
            for (int i = 0; i < result.samples.Length; ++i) result.samples[i] *= gain;
        }
        Renderers.ApplyDynamics(phrase, result);
        progress.Complete(phrase.phones.Length, info);
        return result;
    }, cancellation.Token);

    static void AddDirects(RenderPhrase phrase, ResamplerItem[] items, float[] phraseSamples) {
        foreach (var item in items.Where(item => item.phone.direct)) {
            double posMs = item.phone.positionMs - item.phone.leadingMs -
                (phrase.positionMs - phrase.leadingMs);
            int destination = (int)(posMs * 44.1);
            using var stream = Wave.OpenFile(item.phone.oto.File);
            var samples = Wave.GetSamples(stream.ToSampleProvider().ToMono(1, 0));
            int offset = (int)(item.phone.oto.Offset * 44.1);
            int cutoff = (int)(item.phone.oto.Cutoff * 44.1);
            int length = cutoff >= 0 ? samples.Length - offset - cutoff : -cutoff;
            samples = samples.Skip(offset).Take(length).ToArray();
            item.ApplyEnvelope(samples);
            for (int i = 0; i < samples.Length && destination + i < phraseSamples.Length; ++i) {
                if (destination + i >= 0) phraseSamples[destination + i] = samples[i];
            }
        }
    }
    public RenderPitchResult LoadRenderedPitch(RenderPhrase phrase) => null!;
    // Static suggestions are supplied automatically from plugin metadata by the host.
    public UExpressionDescriptor[] GetSuggestedExpressions(USinger singer, URenderSettings settings) =>
        Array.Empty<UExpressionDescriptor>();
    public override string ToString() => "MORESAMPLER2";
}

sealed class Moresampler2Engine : IResampler {
    static readonly object nativeLock = new();
    [StructLayout(LayoutKind.Sequential)] struct NativeSynthesis {
        public IntPtr input, output; public float tone, velocity; public IntPtr flags;
        public float offset, length, consonant, cutoff; public int volume, modulation, tempo;
        public IntPtr pitchCurve, pitchValues; public int pitchValueCount;
        public IntPtr genderCurve, tensionCurve, breathinessCurve, voicingCurve;
        public IntPtr mouthCurve, formantCurve, singerFormantCurve, growlCurve;
        public IntPtr growlCoarsenessCurve, distortionCurve, pitchRandomnessCurve;
        public IntPtr unvoicedEmphasisCurve, amplitudeModulationCurve;
        public IntPtr stabilityCurve, drynessCurve;
        public int expressionCount;
        public IntPtr analysisFile;
    }
    [StructLayout(LayoutKind.Sequential)] struct NativePhone {
        public NativeSynthesis synthesis; public float positionMs, skipMs, lengthMs, overlapMs;
        public IntPtr envelopeXMs, envelopeY; public int envelopeCount;
    }
    [StructLayout(LayoutKind.Sequential)] struct NativePhrase {
        public IntPtr phones; public int phoneCount; public IntPtr output; public float durationMs; public int sampleRate;
    }
    [StructLayout(LayoutKind.Sequential)] struct NativeResult {
        public IntPtr samples; public int sampleCount, sampleRate;
    }
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)] delegate int RenderPhraseNative(ref NativePhrase request, out NativeResult result);
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)] delegate void FreeNative(IntPtr pointer);
    readonly RenderPhraseNative renderPhrase;
    readonly FreeNative free;
    readonly RendererCacheService cache;
    readonly RendererAnalysisService analysis;
    public string FilePath { get; }
    public bool NoWrapperScript => true;
    public ResamplerManifest Manifest { get; }

    public Moresampler2Engine(RendererPluginContext context) {
        Manifest = context.Manifest;
        cache = context.Cache;
        analysis = context.Analysis;
        var name = OperatingSystem.IsWindows() ? "moresampler2.dll" : OperatingSystem.IsMacOS() ? "libmoresampler2.dylib" : "libmoresampler2.so";
        FilePath = ExtractNativeLibrary(name);
        var library = NativeLibrary.Load(FilePath);
        RendererProfile.Configure(library);
        renderPhrase = Marshal.GetDelegateForFunctionPointer<RenderPhraseNative>(NativeLibrary.GetExport(library, "ms2_render_phrase_v1"));
        free = Marshal.GetDelegateForFunctionPointer<FreeNative>(NativeLibrary.GetExport(library, "ms2_free"));
    }

    public float[] RenderPhrase(RenderPhrase phrase, ResamplerItem[] items, CancellationToken cancellation) {
        using var timing = RendererProfile.Measure("phrase", $"hash={phrase.hash:x16} phones={items.Length}");
        // Stable across renders and rebuilds. This revision excludes pre-fix audio;
        // the host phrase hash distinguishes note/curve/setting changes.
        var outputPath = cache.GetPhrasePath(phrase, ".neutral1.wav");
        if (!phrase.phones.Any(phone => phone.direct) && File.Exists(outputPath)) {
            RendererProfile.Write("phrase.cache_hit", outputPath);
            using var stream = Wave.OpenFile(outputPath);
            return Wave.GetSamples(stream.ToSampleProvider().ToMono(1, 0));
        }
        RendererProfile.Write("phrase.cache_miss_or_direct_bypass", outputPath);
        var preparationTimer = System.Diagnostics.Stopwatch.StartNew();
        using var memory = new NativeMemory();
        int size = Marshal.SizeOf<NativePhone>();
        var phones = Marshal.AllocHGlobal(size * items.Length); memory.Add(phones);
        for (int i = 0; i < items.Length; ++i) {
            cancellation.ThrowIfCancellationRequested();
            var item = items[i];
            var envelope = item.EnvelopeMsToSamples();
            var pitch = item.pitches.Select(value => (float)value).ToArray();
            var native = new NativePhone {
                synthesis = new NativeSynthesis {
                    input = memory.String(item.inputFile), output = IntPtr.Zero,
                    analysisFile = memory.String(analysis.GetPath("llsm2", item.inputFile)),
                    tone = (float)MusicMath.ToneToFreq(item.tone), velocity = item.velocity,
                    flags = memory.String(item.GetFlagsString()), offset = (float)item.offset,
                    length = (float)item.durRequired, consonant = (float)item.consonant, cutoff = (float)item.cutoff,
                    volume = item.volume, modulation = item.modulation, tempo = (int)item.tempo,
                    pitchValues = memory.Floats(pitch), pitchValueCount = pitch.Length,
                    genderCurve = memory.Floats(SampleCurve(phrase, item, phrase.gender)),
                    tensionCurve = memory.Floats(SampleCurve(phrase, item, phrase.tension)),
                    breathinessCurve = memory.Floats(SampleCurve(phrase, item, phrase.breathiness)),
                    voicingCurve = memory.Floats(SampleCurve(phrase, item, phrase.voicing, 100)),
                    mouthCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MOPN"))),
                    formantCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MFRM"))),
                    singerFormantCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MSFM"))),
                    growlCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MGRW"))),
                    growlCoarsenessCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MCRS"))),
                    distortionCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MDST"))),
                    pitchRandomnessCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MPRN"))),
                    unvoicedEmphasisCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MUVC"))),
                    amplitudeModulationCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MAMP"))),
                    stabilityCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MSTA"))),
                    drynessCurve = memory.Floats(SampleCurve(phrase, item, FindCurve(phrase, "MDRY"))),
                    expressionCount = pitch.Length,
                },
                positionMs = (float)(item.phone.positionMs - item.phone.leadingMs - (phrase.positionMs - phrase.leadingMs)),
                skipMs = (float)item.skipOver,
                lengthMs = item.phone.envelope[4].X - item.phone.envelope[0].X,
                overlapMs = item.overlap,
                envelopeXMs = memory.Floats(envelope.Select(p => p.X * 1000f / 44100f).ToArray()),
                envelopeY = memory.Floats(envelope.Select(p => p.Y).ToArray()), envelopeCount = envelope.Count,
            };
            Marshal.StructureToPtr(native, phones + i * size, false);
        }
        var request = new NativePhrase { phones = phones, phoneCount = items.Length,
            output = memory.String(outputPath),
            durationMs = (float)(phrase.durationMs + phrase.leadingMs), sampleRate = 44100 };
        int code; NativeResult result;
        RendererProfile.Write("phrase.prepare", $"elapsed_ms={preparationTimer.Elapsed.TotalMilliseconds:F3} path={outputPath}");
        var lockTimer = System.Diagnostics.Stopwatch.StartNew();
        lock (nativeLock) {
            RendererProfile.Write("phrase.native_lock_wait", $"elapsed_ms={lockTimer.Elapsed.TotalMilliseconds:F3} path={outputPath}");
            using var nativeTiming = RendererProfile.Measure("phrase.native", outputPath);
            code = renderPhrase(ref request, out result);
            RendererProfile.Write("phrase.native_result", $"code={code} samples={result.sampleCount} path={outputPath}");
        }
        if (code != 0 || result.samples == IntPtr.Zero || result.sampleCount <= 0)
            throw new ResamplerFailedException($"moresampler2 phrase renderer returned {code}.");
        try { var samples = new float[result.sampleCount]; Marshal.Copy(result.samples, samples, 0, samples.Length); return samples; }
        finally { free(result.samples); }
    }

    static float[] SampleCurve(RenderPhrase phrase, ResamplerItem item, float[]? curve, float neutral = 0) {
        // Null pointers let the native engine use its neutral defaults directly.
        if (curve == null || curve.Length == 0 || curve.All(value => value == neutral))
            return Array.Empty<float>();
        var result = new float[item.pitches.Length];
        var phraseStart = (int)Math.Floor(phrase.timeAxis.MsPosToNonExactTickPos(phrase.positionMs - phrase.leadingMs));
        var leading = item.phone.oto.Preutter * Math.Pow(2, 1.0 - item.velocity * .01);
        var startMs = item.phone.positionMs - leading;
        var intervalMs = MusicMath.TempoTickToMs(item.tempo, 5);
        for (int i = 0; i < result.Length; ++i) {
            var tick = phrase.timeAxis.MsPosToNonExactTickPos(startMs + intervalMs * i);
            var index = Math.Clamp((tick - phraseStart) / 5.0, 0, curve.Length - 1);
            int left = (int)Math.Floor(index), right = (int)Math.Ceiling(index);
            result[i] = (float)(curve[left] + (curve[right] - curve[left]) * (index - left));
        }
        return result;
    }

    static float[]? FindCurve(RenderPhrase phrase, string abbr) =>
        phrase.curves.FirstOrDefault(curve =>
            string.Equals(curve.Item1, abbr, StringComparison.OrdinalIgnoreCase))?.Item2;

    internal static string ExtractNativeLibrary(string name) {
        var assembly = typeof(Moresampler2Engine).Assembly;
        var resource = assembly.GetManifestResourceNames().FirstOrDefault(x => x.EndsWith(name, StringComparison.OrdinalIgnoreCase))
            ?? throw new DllNotFoundException($"Embedded moresampler2 engine {name} is missing.");
        var directory = Path.Combine(Path.GetTempPath(), "OpenUtau", "moresampler2", assembly.ManifestModule.ModuleVersionId.ToString("N"));
        Directory.CreateDirectory(directory);
        var path = Path.Combine(directory, name);
        if (!File.Exists(path)) { using var input = assembly.GetManifestResourceStream(resource)!; using var output = File.Create(path); input.CopyTo(output); }
        return path;
    }
    public string DoResamplerReturnsFile(ResamplerItem item, ILogger logger) => throw new NotSupportedException("Use phrase rendering.");
    public float[] DoResampler(ResamplerItem item, ILogger logger) => throw new NotSupportedException("Use phrase rendering.");
    public void CheckPermissions() { }
    public bool SupportsFlag(string abbr) => true;
    public override string ToString() => "moresampler2-phrase-v1";

    sealed class NativeMemory : IDisposable {
        readonly List<IntPtr> pointers = new();
        public void Add(IntPtr p) => pointers.Add(p);
        public IntPtr String(string value) {
            var bytes = System.Text.Encoding.UTF8.GetBytes((value ?? "") + '\0');
            var p = Marshal.AllocHGlobal(bytes.Length);
            Marshal.Copy(bytes, 0, p, bytes.Length);
            Add(p);
            return p;
        }
        public IntPtr Floats(float[] values) { if (values.Length == 0) return IntPtr.Zero; var p = Marshal.AllocHGlobal(values.Length * 4); Marshal.Copy(values, 0, p, values.Length); Add(p); return p; }
        public void Dispose() { foreach (var p in pointers) Marshal.FreeHGlobal(p); }
    }
}
