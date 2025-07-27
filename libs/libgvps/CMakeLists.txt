add_library(gvps STATIC
    gvps_full.c
    gvps_obsrv.c
    gvps_sampled.c
    gvps_variable.c
)

target_include_directories(gvps PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})