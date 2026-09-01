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

[ExternalRenderer("io.github.plus-minus-null.moresampler2", "moresampler2")]
public sealed class Moresampler2RendererPlugin : IOpenUtauRendererPlugin {
    static UExpressionDescriptor Curve(string name, string abbr) =>
        new(name, abbr, -100, 100, 0) { type = UExpressionType.Curve };
    public int ApiVersion => ExternalRendererRegistry.ApiVersion;
    public RendererPluginMetadata Metadata { get; } = new() {
        Capabilities = new RendererCapabilitiesManifest { cancellation = false, parallelism = 1 },
        AnalysisFormats = new Dictionary<string, AnalysisFormatManifest> {
            ["llsm2"] = new() { name = "moresampler2 LLSM analysis", path = "{wav_dir}/{wav_stem}.llsm2", required = true, canGenerate = true },
        },
        Expressions = new Dictionary<string, UExpressionDescriptor> {
            [Ustx.GENC] = Curve("Gender", Ustx.GENC),
            [Ustx.TENC] = Curve("Tension", Ustx.TENC),
        },
    };
    public IRenderer CreateRenderer(RendererPluginContext context) => new Moresampler2Renderer(context);
    public IRendererAnalysisProvider CreateAnalysisProvider(RendererPluginContext context) =>
        new Moresampler2AnalysisProvider(context.Analysis);
}

sealed class Moresampler2AnalysisProvider : IRendererAnalysisProvider, IDisposable {
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate int AnalyzeNative(IntPtr input, IntPtr output);
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate int ValidateNative(IntPtr analysis);
    readonly IntPtr library;
    readonly AnalyzeNative analyze;
    readonly ValidateNative validate;
    readonly RendererAnalysisService analysis;

    public Moresampler2AnalysisProvider(RendererAnalysisService analysis) {
        this.analysis = analysis;
        var name = OperatingSystem.IsWindows() ? "moresampler2.dll" :
            OperatingSystem.IsMacOS() ? "libmoresampler2.dylib" : "libmoresampler2.so";
        library = NativeLibrary.Load(Moresampler2Engine.ExtractNativeLibrary(name));
        analyze = Marshal.GetDelegateForFunctionPointer<AnalyzeNative>(
            NativeLibrary.GetExport(library, "ms2_analyze_v1"));
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
                cancellation.ThrowIfCancellationRequested();
                var input = Marshal.StringToCoTaskMemUTF8(request.SourceFile);
                var output = Marshal.StringToCoTaskMemUTF8(request.OutputFile);
                try {
                    int code = analyze(input, output);
                    results.Add(code == 0
                        ? new RendererAnalysisResult(request, RendererAnalysisOutcome.Generated)
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
        cancellation.ThrowIfCancellationRequested();
        var basic = analysis.GetBasicState(request.Format, request.SourceFile);
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
        d.abbr is Ustx.DYN or Ustx.PITD or Ustx.VEL or Ustx.VOL or Ustx.MOD or Ustx.GENC or Ustx.TENC;
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
        Renderers.ApplyDynamics(phrase, result);
        progress.Complete(phrase.phones.Length, info);
        PlaybackManager.Inst.LiveWaveformCache[phrase.hash.ToString()] =
            (trackNo, phrase.positionMs - phrase.leadingMs, result.samples, DateTime.Now);
        DocManager.Inst.ExecuteCmd(new WaveformReadyNotification());
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
    public override string ToString() => "moresampler2";
}

sealed class Moresampler2Engine : IResampler {
    static readonly object nativeLock = new();
    [StructLayout(LayoutKind.Sequential)] struct NativeSynthesis {
        public IntPtr input, output; public float tone, velocity; public IntPtr flags;
        public float offset, length, consonant, cutoff; public int volume, modulation, tempo;
        public IntPtr pitchCurve, pitchValues; public int pitchValueCount;
        public IntPtr genderCurve, tensionCurve; public int expressionCount;
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
    public string FilePath { get; }
    public ResamplerManifest Manifest { get; }

    public Moresampler2Engine(RendererPluginContext context) {
        Manifest = context.Manifest;
        cache = context.Cache;
        var name = OperatingSystem.IsWindows() ? "moresampler2.dll" : OperatingSystem.IsMacOS() ? "libmoresampler2.dylib" : "libmoresampler2.so";
        FilePath = ExtractNativeLibrary(name);
        var library = NativeLibrary.Load(FilePath);
        renderPhrase = Marshal.GetDelegateForFunctionPointer<RenderPhraseNative>(NativeLibrary.GetExport(library, "ms2_render_phrase_v1"));
        free = Marshal.GetDelegateForFunctionPointer<FreeNative>(NativeLibrary.GetExport(library, "ms2_free"));
    }

    public float[] RenderPhrase(RenderPhrase phrase, ResamplerItem[] items, CancellationToken cancellation) {
        var outputPath = cache.GetPhrasePath(phrase);
        if (!phrase.phones.Any(phone => phone.direct) && File.Exists(outputPath)) {
            using var stream = Wave.OpenFile(outputPath);
            return Wave.GetSamples(stream.ToSampleProvider().ToMono(1, 0));
        }
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
                    tone = (float)MusicMath.ToneToFreq(item.tone), velocity = item.velocity,
                    flags = memory.String(item.GetFlagsString()), offset = (float)item.offset,
                    length = (float)item.durRequired, consonant = (float)item.consonant, cutoff = (float)item.cutoff,
                    volume = item.volume, modulation = item.modulation, tempo = (int)item.tempo,
                    pitchValues = memory.Floats(pitch), pitchValueCount = pitch.Length,
                    genderCurve = memory.Floats(SampleCurve(phrase, item, phrase.gender)),
                    tensionCurve = memory.Floats(SampleCurve(phrase, item, phrase.tension)), expressionCount = pitch.Length,
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
        lock (nativeLock) code = renderPhrase(ref request, out result);
        if (code != 0 || result.samples == IntPtr.Zero || result.sampleCount <= 0)
            throw new ResamplerFailedException($"moresampler2 phrase renderer returned {code}.");
        try { var samples = new float[result.sampleCount]; Marshal.Copy(result.samples, samples, 0, samples.Length); return samples; }
        finally { free(result.samples); }
    }

    static float[] SampleCurve(RenderPhrase phrase, ResamplerItem item, float[]? curve) {
        var result = new float[item.pitches.Length];
        if (curve == null || curve.Length == 0) return result;
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
