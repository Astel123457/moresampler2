# moresampler2
moresampler2, written in c, with libllsm2 support

For regular OpenUtau CLASSIC use, place [`moresampler2.yaml`](moresampler2.yaml)
beside the standalone `moresampler2.exe` (or `moresampler2` on Linux/macOS) in
Resamplers, restart OpenUtau, and select it in CLASSIC's track settings.
This manifest exposes per-note flags without requiring an external-renderer API
build or plugin DLL. If you rename the executable, rename the YAML to match its
basename. It does not add renderer-only curves or analysis-generation controls.

An optional managed bridge for OpenUtau's external-renderer discovery API is in
[`OpenUtau/`](OpenUtau/). It embeds the native engine in one platform-specific
DLL while preserving the standalone executable for normal resampler compatibility.

Listen, I don't know c. If the code looks horrible, it prolly is, i just changed things till they worked, with a mix of chatgpt and screaming

I can't guarntee that everthing works on mac and linux, i use windows almost exlusively

libllsm2 is finnicky, and also very slow at analysis unfortunately, until openutau supports the llsm2 cache files it will take a while to render,
i suggest just running moresampler2 in the background for a bit before rendering, as the actual render times of the files is suprisingly fast

also moresampler2 is basically in beta, don't expect everything from original moresampler to work directly in moresampler 2 until it's implemented

moresampler2 currently can't
1. auto oto
2. consonant velocity works but is broken
3. no moreconfig.txt to change options on

i understand if you want certain things now, however. i'm only a guy, and unless you can code it faster than me, you'll just have to wait a little while in-between updates
