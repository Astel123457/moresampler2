cmake_minimum_required(VERSION 3.10)

file(READ "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION_FROM_FILE)
string(STRIP "${PROJECT_VERSION_FROM_FILE}" PROJECT_VERSION_FROM_FILE)

project(moresampler2 VERSION ${PROJECT_VERSION_FROM_FILE} LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_definitions(-DFP_TYPE=float)

add_subdirectory(libs/libgvps)
add_subdirectory(libs/libpyin)
add_subdirectory(libs/ciglet)
add_subdirectory(libs/libllsm)

add_executable(moresampler2 moresampler2.c wavtool.c utils.c)

add_library(moresampler2_engine SHARED moresampler2.c wavtool.c utils.c)
set_target_properties(moresampler2_engine PROPERTIES OUTPUT_NAME moresampler2)
target_compile_definitions(moresampler2_engine PRIVATE
    MORESAMPLER2_NO_MAIN
    MORESAMPLER2_BUILD_SHARED
    MORESAMPLER2_VERSION="${PROJECT_VERSION_FROM_FILE}"
)

target_compile_definitions(moresampler2 PRIVATE
    MORESAMPLER2_VERSION="${PROJECT_VERSION_FROM_FILE}"
)

target_link_libraries(moresampler2 PRIVATE
    llsm
    ciglet
    pyin
    gvps
)

target_link_libraries(moresampler2_engine PRIVATE
    llsm
    ciglet
    pyin
    gvps
)

if(NOT WIN32)
    target_link_libraries(moresampler2 PRIVATE m)
    target_link_libraries(moresampler2_engine PRIVATE m)
endif()

target_include_directories(moresampler2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libllsm
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/ciglet
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpyin
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libgvps
)

target_include_directories(moresampler2_engine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libllsm
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/ciglet
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpyin
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libgvps
)

include(CTest)
option(MS2_CHECK_ALLOCATIONS "Check native test allocation balance (Linux/GNU linker)" OFF)
if(BUILD_TESTING)
    add_executable(test_velocity tests/velocity.c utils.c wavtool.c)
    target_link_libraries(test_velocity PRIVATE llsm ciglet pyin gvps)
    if(NOT WIN32)
        target_link_libraries(test_velocity PRIVATE m)
    endif()
    add_test(NAME velocity_mapping COMMAND test_velocity)
    add_executable(test_stitch tests/stitch.c utils.c)
    target_link_libraries(test_stitch PRIVATE llsm ciglet pyin gvps)
    if(NOT WIN32)
        target_link_libraries(test_stitch PRIVATE m)
    endif()
    add_test(NAME phrase_stitch COMMAND test_stitch)
    add_executable(test_defaults tests/defaults.c utils.c wavtool.c)
    target_link_libraries(test_defaults PRIVATE llsm ciglet pyin gvps)
    if(NOT WIN32)
        target_link_libraries(test_defaults PRIVATE m)
    endif()
    add_test(NAME neutral_effect_defaults COMMAND test_defaults)
    add_executable(test_phrase_cache tests/phrase_cache.c utils.c wavtool.c)
    target_link_libraries(test_phrase_cache PRIVATE llsm ciglet pyin gvps)
    if(NOT WIN32)
        target_link_libraries(test_phrase_cache PRIVATE m)
    endif()
    add_test(NAME phrase_cache_pipeline COMMAND test_phrase_cache)
    if(MS2_CHECK_ALLOCATIONS)
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
            message(FATAL_ERROR "MS2_CHECK_ALLOCATIONS requires Linux and a GNU-compatible linker")
        endif()
        foreach(test_target test_velocity test_stitch test_defaults test_phrase_cache)
            target_sources(${test_target} PRIVATE tests/allocations.c)
            set_property(TARGET ${test_target} APPEND_STRING PROPERTY LINK_FLAGS
                " -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc,--wrap=free")
        endforeach()
    endif()
endif()
