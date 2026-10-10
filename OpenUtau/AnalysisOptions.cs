using System.Globalization;
using System.Runtime.InteropServices;

namespace Moresampler2.OpenUtau;

[StructLayout(LayoutKind.Sequential)]
internal struct AnalysisOptions {
    public float MinimumF0, MaximumF0, VoicingBias;
    public int TransitionRange, RefineF0;

    public static AnalysisOptions From(IReadOnlyDictionary<string, string> settings) {
        float Number(string key, float fallback, float min, float max) {
            if (!settings.TryGetValue(key, out var text)) return fallback;
            if (!float.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out var value)
                    || !float.IsFinite(value) || value < min || value > max)
                throw new InvalidDataException($"Invalid moresampler2 setting {key}: expected {min}–{max}.");
            return value;
        }
        var result = new AnalysisOptions {
            MinimumF0 = Number("analysisF0Min", 50, 20, 1999),
            MaximumF0 = Number("analysisF0Max", 800, 21, 2000),
            VoicingBias = Number("analysisVoicingBias", 2, .1f, 10),
            TransitionRange = (int)Number("analysisTransitionRange", 24, 1, 48),
            RefineF0 = 1,
        };
        if (Number("analysisTransitionRange", 24, 1, 48) != result.TransitionRange)
            throw new InvalidDataException("moresampler2 pitch transition range must be an integer.");
        if (result.MinimumF0 >= result.MaximumF0)
            throw new InvalidDataException("moresampler2 minimum F0 must be lower than maximum F0.");
        if (settings.TryGetValue("analysisRefineF0", out var refine)) {
            if (!bool.TryParse(refine, out var enabled))
                throw new InvalidDataException("moresampler2 F0 refinement must be true or false.");
            result.RefineF0 = enabled ? 1 : 0;
        }
        return result;
    }

}
