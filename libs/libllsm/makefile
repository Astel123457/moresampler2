PREFIX = /usr
CC ?= gcc
AR ?= ar

OUT_DIR = ./build
OBJS = $(OUT_DIR)/container.o \
CORE_OBJS = $(OUT_DIR)/container.o \
            $(OUT_DIR)/frame.o \
            $(OUT_DIR)/dsputils.o \
            $(OUT_DIR)/llsmutils.o \
            $(OUT_DIR)/layer0.o \
            $(OUT_DIR)/layer1.o \
            $(OUT_DIR)/coder.o \
            $(OUT_DIR)/llsmrt.o
TARGET_A = $(OUT_DIR)/libllsm2.a

CIGLET_PREFIX = ./ext
GVPS_PREFIX = ./ext
PYIN_PREFIX = ./ext

FP_TYPE ?= float
CONFIG  ?= Debug

CIGLET_A = $(CIGLET_PREFIX)/libciglet.a
CIGLET_INCLUDE = $(CIGLET_PREFIX)/
GVPS_A = $(GVPS_PREFIX)/libgvps.a
GVPS_INCLUDE = $(GVPS_PREFIX)/
PYIN_A = $(PYIN_PREFIX)/libpyin.a
PYIN_INCLUDE = $(PYIN_PREFIX)/

ARFLAGS = -rv
CFLAGS_COMMON = -DFP_TYPE=$(FP_TYPE) -std=c99 -Wall -I$(CIGLET_INCLUDE) -I$(GVPS_INCLUDE) -I$(PYIN_INCLUDE)
ifeq ($(CXX), emcc)
  CFLAGS_DBG = $(CFLAGS_COMMON) -O1 -g -D_DEBUG
  CFLAGS_REL = $(CFLAGS_COMMON) -O3
else
  CFLAGS_COMMON += -fPIC -pthread -DUSE_PTHREAD
  CFLAGS_DBG = $(CFLAGS_COMMON) -fopenmp -Og -g -D_DEBUG
  CFLAGS_REL = $(CFLAGS_COMMON) -fopenmp -Ofast
endif
ifeq ($(CONFIG), Debug)
  CFLAGS = $(CFLAGS_DBG)
else
  CFLAGS = $(CFLAGS_REL)
endif

default: $(TARGET_A)

test: test-dsputils test-layer0 test-layer1 test-coder test-llsmrt test-pbpeffects

ifeq ($(CXX), emcc)

test-dsputils: $(OUT_DIR)/test-structs.js \
               $(OUT_DIR)/test-dsputils.js \
               $(OUT_DIR)/test-harmonic.js
	node --experimental-wasm-threads $(OUT_DIR)/test-structs.js
	node --experimental-wasm-threads $(OUT_DIR)/test-dsputils.js
	node --experimental-wasm-threads $(OUT_DIR)/test-harmonic.js

test-layer0: $(OUT_DIR)/test-layer0-anasynth.html \
             $(OUT_DIR)/test-layer0-edgecase.html
	emrun $(OUT_DIR)/test-layer0-anasynth.html pp
	emrun $(OUT_DIR)/test-layer0-anasynth.html czt
	emrun $(OUT_DIR)/test-layer0-edgecase.html

test-layer1: $(OUT_DIR)/test-layer1-anasynth.html
	emrun $(OUT_DIR)/test-layer1-anasynth.html

test-coder: $(OUT_DIR)/test-coder.html
	emrun $(OUT_DIR)/test-coder.html

test-llsmrt: $(OUT_DIR)/test-llsmrt.html
	emrun $(OUT_DIR)/test-llsmrt.html

test-pbpeffects: $(OUT_DIR)/test-pbpeffects.html
	emrun $(OUT_DIR)/test-pbpeffects.html

else

test-dsputils: $(OUT_DIR)/test-structs.exe \
               $(OUT_DIR)/test-dsputils.exe \
               $(OUT_DIR)/test-harmonic.exe
	$(OUT_DIR)/test-structs.exe
	$(OUT_DIR)/test-dsputils.exe
	$(OUT_DIR)/test-harmonic.exe

test-layer0: $(OUT_DIR)/test-layer0-anasynth.exe \
             $(OUT_DIR)/test-layer0-edgecase.exe
	$(OUT_DIR)/test-layer0-anasynth.exe pp
	$(OUT_DIR)/test-layer0-anasynth.exe czt
	$(OUT_DIR)/test-layer0-edgecase.exe

test-layer1: $(OUT_DIR)/test-layer1-anasynth.exe
	$(OUT_DIR)/test-layer1-anasynth.exe

test-coder: $(OUT_DIR)/test-coder.exe
	$(OUT_DIR)/test-coder.exe

test-llsmrt: $(OUT_DIR)/test-llsmrt.exe
	$(OUT_DIR)/test-llsmrt.exe

test-pbpeffects: $(OUT_DIR)/test-pbpeffects.exe
	$(OUT_DIR)/test-pbpeffects.exe

endif


$(OUT_DIR)/test-%: test/test-%.c $(TARGET_A) \
	  $(CIGLET_A) $(PYIN_A) $(GVPS_A) test/verify-utils.h
	$(CC) $(CFLAGS) -o $(OUT_DIR)/test-$* test/test-$*.c \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(GVPS_A) -lm

$(OUT_DIR)/test-%.js: test/test-%.c $(TARGET_A) \
	  $(CIGLET_A) $(PYIN_A) $(GVPS_A) test/verify-utils.h
	$(CC) $(CFLAGS) -o $(OUT_DIR)/test-$*.js test/test-$*.c \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(GVPS_A) -lm

$(OUT_DIR)/test-%.html: test/test-%.c $(TARGET_A) \
	  $(CIGLET_A) $(PYIN_A) $(GVPS_A) test/verify-utils.h
	$(CC) $(CFLAGS) -o $(OUT_DIR)/test-$*.html test/test-$*.c \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(GVPS_A) -lm \
	  --preload-file test/arctic_a0001.wav \
	  --preload-file test/are-you-ready.wav --emrun \
	  -s TOTAL_MEMORY=128MB -s PROXY_TO_PTHREAD

$(OUT_DIR)/demo-%: test/demo-%.c $(TARGET_A) \
	  $(CIGLET_A) $(PYIN_A) $(GVPS_A)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/demo-$* test/demo-$*.c $(CFLAGS_PLAT) \
	  $(TARGET_A) $(CIGLET_A) $(PYIN_A) $(GVPS_A) -lm

$(TARGET_A): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET_A) $(OBJS)

$(OUT_DIR)/frame.o: llsm.h dsputils.h
$(OUT_DIR)/container.o: llsm.h
$(OUT_DIR)/dsputils.o: dsputils.h llsm.h
$(OUT_DIR)/llsmutils.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/layer0.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/layer1.o: llsmutils.h dsputils.h llsm.h
$(OUT_DIR)/coder.o: dsputils.h llsm.h
$(OUT_DIR)/llsmrt.o: buffer.h llsmutils.h dsputils.h llsm.h llsmrt.h

$(OUT_DIR)/%.o: %.c constants.h
	mkdir -p $(OUT_DIR)
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c $*.c

install: $(OUT_DIR)/libllsm2.a 
	mkdir -p $(PREFIX)/lib $(PREFIX)/include/libllsm2
	cp $(OUT_DIR)/libllsm2.a $(PREFIX)/lib
	cp llsm.h llsmrt.h llsmutils.h dsputils.h buffer.h $(PREFIX)/include/libllsm2

clean:
	rm -f build/*

$(OUT_DIR)/test-structs.exe: test/test-structs.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-dsputils.exe: test/test-dsputils.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-harmonic.exe: test/test-harmonic.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-layer0-anasynth.exe: test/test-layer0-anasynth.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-layer0-edgecase.exe: test/test-layer0-edgecase.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-layer1-anasynth.exe: test/test-layer1-anasynth.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-coder.exe: test/test-coder.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-llsmrt.exe: test/test-llsmrt.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)

$(OUT_DIR)/test-pbpeffects.exe: test/test-pbpeffects.c
	$(CC) $(CFLAGS) $< -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)


moresampler2: $(OUT_DIR)/moresampler2.exe

$(OUT_DIR)/moresampler2.exe: moresampler2.c $(CORE_OBJS)
	$(CC) $(CFLAGS) moresampler2.c $(CORE_OBJS) -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)
	$(CC) $(CFLAGS) $^ -o $@ $(CIGLET_A) $(GVPS_A) $(PYIN_A)