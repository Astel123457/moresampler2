add_library(pyin STATIC
    pyin.c
    yin.c
    math-funcs.c
    matlabfunctions.c
)

target_link_libraries(pyin PUBLIC gvps)

target_include_directories(pyin PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../libgvps
)