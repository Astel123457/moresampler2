# OpenUtau renderer bridge

This bridge exposes moresampler2 through OpenUtau's external-renderer discovery
API. `Moresampler2Renderer` receives complete OpenUtau phrases and calls the
native engine in-process for synthesis. The original executable remains available
for normal UTAU resampler compatibility.

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

The plugin class declares its identity, expressions, capabilities, and `.llsm2`
analysis relationship. No YAML file is required. This produces a platform-specific
bridge DLL because its matching native engine is embedded inside it.
