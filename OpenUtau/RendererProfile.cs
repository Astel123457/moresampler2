using System.Diagnostics;
using System.Runtime.InteropServices;
using OpenUtau.Core;

namespace Moresampler2.OpenUtau;

// Temporary diagnostic profiling; one pair of append-only logs per app process.
internal static class RendererProfile {
    static readonly object gate = new();
    static readonly HashSet<IntPtr> configured = new();
    static readonly string prefix = Path.Combine(PathManager.Inst.LogsPath,
        $"moresampler2-profile-{DateTime.UtcNow:yyyyMMdd-HHmmss}-{Environment.ProcessId}");
    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    delegate void SetPath(IntPtr path);

    public static void Configure(IntPtr library) {
        lock (gate) {
            if (configured.Contains(library)) return;
            try {
                Directory.CreateDirectory(PathManager.Inst.LogsPath);
                var path = Marshal.StringToCoTaskMemUTF8(prefix + ".native.log");
                try {
                    Marshal.GetDelegateForFunctionPointer<SetPath>(
                        NativeLibrary.GetExport(library, "ms2_set_profile_path_v1"))(path);
                } finally { Marshal.FreeCoTaskMem(path); }
                configured.Add(library);
                Write("session", $"assembly={typeof(RendererProfile).Assembly.ManifestModule.ModuleVersionId} native={prefix}.native.log");
            } catch (Exception e) {
                Serilog.Log.Warning(e, "Could not enable moresampler2 profiling");
            }
        }
    }

    public static void Write(string stage, string detail) {
        try {
            lock (gate) File.AppendAllText(prefix + ".managed.log",
                $"{DateTime.UtcNow:O} thread={Environment.CurrentManagedThreadId} stage={stage} {detail.Replace('\r', ' ').Replace('\n', ' ')}{Environment.NewLine}");
        } catch { /* Diagnostics must not fail synthesis. */ }
    }

    public static IDisposable Measure(string stage, string detail) => new Timing(stage, detail);
    sealed class Timing : IDisposable {
        readonly string stage, detail;
        readonly Stopwatch timer = Stopwatch.StartNew();
        public Timing(string stage, string detail) {
            this.stage = stage; this.detail = detail;
            Write(stage + ".begin", detail);
        }
        public void Dispose() => Write(stage + ".end",
            $"elapsed_ms={timer.Elapsed.TotalMilliseconds:F3} {detail}");
    }
}
