add_library(ciglet STATIC
	ciglet.c
	fast_median.c
	wavfile.c
	fftsg_h.c
)

target_include_directories(ciglet PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})