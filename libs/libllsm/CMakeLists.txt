add_library(llsm STATIC
    container.c
    frame.c
    dsputils.c
    llsmutils.c
    layer0.c
    layer1.c
    coder.c
    llsmrt.c
)

target_link_libraries(llsm PUBLIC ciglet)
target_include_directories(llsm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../
)