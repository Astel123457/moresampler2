# OpenUtau renderer bridge

This bridge exposes moresampler2 through OpenUtau's external-renderer discovery
API. `Moresampler2Renderer` receives complete OpenUtau phrases and calls the
native engine in-process for synthesis. The original executable remains available
for normal UTAU resampler compatibility.

Optional effects are bypassed at neutral settings. Missing/flat curves do not
enable effects; voicing is neutral at 100, and noise suppression (`MSY` / `Y`)
is neutral at 100. Existing projects created with the former MSY default of 0
may retain that descriptor or curve: set it to 100 to restore neutral behavior.
An explicit `Y0` remains an active effect.

Peak normalization is now opt-in (`P0` by default); use `P86` explicitly to
restore the previous automatic level adjustment. Neutral output can therefore
be quieter. Source analysis, time mapping, target pitch and pitch curves,
pitch-related level compensation, volume, envelopes, and stitching still run
as required for resampling. This is not a raw-waveform passthrough.

Renderer mode keeps per-phone model frames in memory and writes only the final
stitched phrase WAV, plus reusable source `.llsm2` analysis when needed. With
`P0`, it does not synthesize per-phone PCM. Explicit `P` normalization may create
temporary PCM in memory to measure peaks, but never writes intermediate WAVs
or resampled `.llsm2` files. The standalone UTAU resampler retains its file-based
interface.

Analysis and phrase caches are reused across renders and rebuilds. Phrase cache
keys follow the host's phrase hash so changed notes, curves, or settings do not
reuse old audio. Automatic rendering generates missing or unusable analysis.
Valid analysis is reused regardless of source timestamps or changed analysis settings.
Unreadable/unsupported analysis is regenerated automatically with a warning in
the logs. Errors are reported only if regeneration fails.

Source analysis is stored as Layer 1, including the glottal and spectral analysis
needed for rendering. The Singers-window generation API performs this conversion
before saving. Old Layer-0 caches are automatically regenerated in renderer mode.
The standalone resampler can convert legacy data in memory without changing the
file. Explicit Singers-window generation overwrites the selected analysis using
default settings. Saves use a temporary file and replacement, so failed writes
do not truncate the existing cache.
Phase propagation remains render-local and is not baked repeatedly into analysis.

### Track settings

Alongside output gain, the renderer exposes:

| Setting | Default | Meaning |
| --- | --- | --- |
| Minimum / maximum detected F0 | 50 / 800 Hz | pYIN source pitch search interval; minimum must be below maximum. |
| Analysis voicing bias | 2 | Above 1 favors voiced frames; below 1 favors unvoiced frames. |
| Pitch transition range | 24 bins | Maximum tracker transition distance in pYIN's quantized pitch bins, not semitones. |
| Refine detected F0 | Enabled | libllsm2 refines the detected source F0. |

Analysis settings affect newly generated files only. All tracks share the source's
existing `.llsm2`, even if their analysis settings differ; the settings used for
the first generation remain effective until the user overwrites or removes it.
Changing settings does not create a separate cache variant. Singers-window batch
generation explicitly overwrites selected files using defaults; to generate with
track settings, remove the selected source analysis and render that track.
These are pYIN/libllsm2 controls, not a
claim of exact equivalence to the original moresampler's analysis algorithm.
Original-only options such as FRQ import and analysis noise reduction are not
exposed until their processing is implemented.

The host passes track settings in `RendererAnalysisRequest.Settings` and does not
authorize overwriting during automatic preparation. Moresampler2 owns the reuse
policy and uses the standard source-analysis path. This bridge must be built with the
matching native library and updated host API (its native synthesis request now
includes an optional analysis-file path).

### Profiling a slow render

This diagnostic build writes a new pair of append-only logs per OpenUtau process
in OpenUtau's `Logs` directory: `moresampler2-profile-<UTC date-time>-<pid>.managed.log`
and `.native.log`. Restart OpenUtau, reproduce the slow render, and keep both files.
Profiling does not clear or rename audio/analysis caches. To measure uncached
synthesis, explicitly rerender the project using OpenUtau's render controls;
do not delete the singer's analysis files.

Managed records include UTC timestamps, phrase hashes, cache hits, validation
states, generation results, preparation time, and native-lock wait time. Native
records include monotonic timestamps and wall-clock milliseconds for analysis
reads/generation, source conversion, time mapping, effects, per-phone synthesis,
accumulated-phrase copying, stitching, final synthesis, and WAV writes. Stage
totals are nested, so do not sum totals and their sub-stages together. Begin
records without matching completion records help locate interrupted work.
Logs include local file paths. Profiling adds small logging overhead and should
be disabled/removed after the performance investigation.

Build it against a checkout containing the external renderer API:

```sh
dotnet build OpenUtau/Moresampler2.OpenUtau.csproj -c Release \
  -p:OpenUtauCoreProject=/path/to/OpenUtau/OpenUtau.Core/OpenUtau.Core.csproj
```

Build the native executable and shared library with CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

The native build must exist before building the managed bridge so it can be
embedded. Copy the resulting single file into OpenUtau's `Resamplers` directory:

- `OpenUtau/bin/Release/net10.0/Moresampler2.OpenUtau.dll`

Do not distribute a private `OpenUtau.Core.dll`. OpenUtau supplies the shared
host assembly when it loads the bridge.

While the external renderer API is in development, rebuild this bridge against
the OpenUtau version you run after updating the host. Even a field changing to
a property (such as `RenderPhrase.phones`) requires recompilation; an older DLL
can otherwise fail with `MissingFieldException`. Restart OpenUtau after replacing
the DLL. Waveform publication is handled by the host from the returned render result.

The plugin class declares its identity, expressions, capabilities, and `.llsm2`
analysis relationship. No YAML file is required. This produces a platform-specific
bridge DLL because its matching native engine is embedded inside it.
