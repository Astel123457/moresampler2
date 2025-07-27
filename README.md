# moresampler2
moresampler2, written in c, with libllsm2 support

Listen, I don't know c. If the code looks horrible, it prolly is, i just changed things till they worked, with a mix of chatgpt and screaming

I can't guarntee that everthing works on mac and linux, i use windows almost exlusively

libllsm2 is finnicky, and also very slow at analysis unfortunately, until openutau supports the llsm2 cache files it will take a while to render,
i suggest just running moresampler2 in the background for a bit before rendering, as the actual render times of the files is suprisingly fast

also moresampler2 is basically in beta, don't expect everything from original moresampler to work directly in moresampler 2 until it's implemented

moresampler2 currently can't
1. auto oto
2. do anything with flags
3. wavtool
4. consonant velocity doesn't do anything
5. no moreconfig.txt to change options on
6. basically anything that sets moresampler apart from the rest (other than using moresamplers second verion of it's engine)

i understand if you want certain things now, however. i'm only a guy, and unless you can code it faster than me, you'll just have to wait a little while in-between updates