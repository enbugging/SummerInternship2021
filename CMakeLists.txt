cmake_minimum_required(VERSION 3.14)                        # CMake version check

project(CHARMM_Optimizer                                    # Set the project name
        VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 11)                                  # Enable C++11 standard. Note GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD_REQUIRED ON)                         # Enforce C++11 standard.

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR                  # Enforce use of GNU Compiler for expected behaviours
    CMAKE_CXX_COMPILER_ID MATCHES "AppleClang")

    message("Build from stable code base")
    file(GLOB_RECURSE 
        SOURCES 
        "lib/**.c"                                          # External librarbies source files
        "lib/**.cpp" 
        "src/**.cpp")                                       # Stable source files
    file(GLOB_RECURSE 
        HEADERS 
        "lib/**.h"                                          # External libraries header files
        "lib/**.hpp" 
        "src/**.hpp")                                       # Stable header files
    add_executable(CHARMM_Optimizer ${SOURCES} ${HEADERS})     # Source files
endif()
