using System.Globalization;
using Moresampler2.OpenUtau;

static void Check(bool ok, string message) {
    if (!ok) throw new Exception(message);
}
var settings = new Dictionary<string, string>();
var defaults = AnalysisOptions.From(settings);
Check(defaults.MinimumF0 == 50 && defaults.MaximumF0 == 800, "Default F0 range changed");
settings["outputGainDb"] = "6";
Check(AnalysisOptions.From(settings).Equals(defaults), "Output gain must not change analysis options");
foreach (var pair in new Dictionary<string, string> {
    ["analysisF0Min"] = "60", ["analysisF0Max"] = "900",
    ["analysisVoicingBias"] = "1.5", ["analysisTransitionRange"] = "12",
    ["analysisRefineF0"] = "false",
}) {
    settings[pair.Key] = pair.Value;
    var configured = AnalysisOptions.From(settings);
    Check(!configured.Equals(defaults), "Each analysis option must be applied to new generation");
    settings.Remove(pair.Key);
}
settings["analysisVoicingBias"] = "1.5";
var before = AnalysisOptions.From(settings);
CultureInfo.CurrentCulture = CultureInfo.GetCultureInfo("de-DE");
Check(AnalysisOptions.From(settings).Equals(before), "Analysis settings must be culture-independent");
foreach (var invalid in new[] {
    new Dictionary<string,string> { ["analysisF0Min"] = "900" },
    new Dictionary<string,string> { ["analysisF0Max"] = "NaN" },
    new Dictionary<string,string> { ["analysisVoicingBias"] = "0" },
    new Dictionary<string,string> { ["analysisTransitionRange"] = "1.5" },
    new Dictionary<string,string> { ["analysisRefineF0"] = "maybe" },
}) {
    bool rejected = false;
    try { AnalysisOptions.From(invalid); } catch (InvalidDataException) { rejected = true; }
    Check(rejected, "Invalid settings must be rejected");
}
Console.WriteLine("Managed analysis-settings checks passed.");
