cmake_minimum_required(VERSION 3.10)
project(moresampler2 C)

add_definitions(-DFP_TYPE=float)

set(CMAKE_C_STANDARD 99)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add all library subdirectories (in dependency order!)
add_subdirectory(libs/libgvps)
add_subdirectory(libs/libpyin)
add_subdirectory(libs/ciglet)
add_subdirectory(libs/libllsm)

# Main program
add_executable(moresampler2 moresampler2.c)

# Link against all the libraries
target_link_libraries(moresampler2 PRIVATE
    llsm
    ciglet
    pyin
    gvps
)

# Let main program see all library headers
target_include_directories(moresampler2 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libllsm
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/ciglet
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libpyin
    ${CMAKE_CURRENT_SOURCE_DIR}/libs/libgvps
)