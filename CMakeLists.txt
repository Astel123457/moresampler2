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

add_executable(moresampler2 moresampler2.c)

target_compile_definitions(moresampler2 PRIVATE
    MORESAMPLER2_VERSION="${PROJECT_VERSION_FROM_FILE}"
)

target_link_libraries(moresampler2 PRIVATE
    llsm
    ciglet
    pyin
    gvps
)

if(NOT WIN32)
    target_link_libraries(moresampler2 PRIVATE m)
endif()

target_include_directories(moresampler2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libllsm
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/ciglet
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpyin
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libgvps
)