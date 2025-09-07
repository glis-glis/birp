cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Policies (for properly setting the timestamps of extracted contents)
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif (POLICY CMP0135)

# To avoid warnings when compaling on MacOS on CRAN
set(CMAKE_OSX_DEPLOYMENT_TARGET "13.7" CACHE STRING "Minimum OS X deployment version")

# Project
project(birp LANGUAGES CXX)

# Needed for Language Server Protocol
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Download coretools (if necessary)
include(FetchContent)
set(CORETOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coretools")

if(EXISTS "${CORETOOLS_DIR}/CMakeLists.txt")
    message(STATUS "Using vendored coretools")
    add_subdirectory(${CORETOOLS_DIR})
else()
    message(STATUS "Fetching coretools")
    FetchContent_Declare(coretools
            URL https://bitbucket.org/wegmannlab/coretools/get/ba8ed2442228a41d65f55cf3007c60bf28da455d.tar.gz
    )
    FetchContent_GetProperties(coretools)
    if(NOT coretools_POPULATED)
        FetchContent_Populate(coretools)
        # Print debugging stuff
        message(STATUS "coretools_SOURCE_DIR = ${coretools_SOURCE_DIR}")
        message(STATUS "coretools_BINARY_DIR = ${coretools_BINARY_DIR}")
        # Copy into libs/ so it gets vendored for CRAN later
        file(COPY ${coretools_SOURCE_DIR}/ DESTINATION ${CORETOOLS_DIR})
    endif()
    add_subdirectory(${CORETOOLS_DIR})
endif()

get_property(_allTargets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "All targets at this point (after coretools): ${_allTargets}")

# Fail now if coretools was not found!
if(NOT TARGET coretools)
    message(FATAL_ERROR "coretools target was not created! Check FetchContent / vendoring setup.")
endif()

# Download stattools (if necessary)
set(STATTOOLS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/stattools")

if(EXISTS "${STATTOOLS_DIR}/CMakeLists.txt")
    message(STATUS "Using vendored stattools")
    add_subdirectory(${STATTOOLS_DIR})
else()
    message(STATUS "Fetching stattools")
    FetchContent_Declare(stattools
            URL https://bitbucket.org/wegmannlab/stattools/get/bae5da94b063f1fce82ebf33b2c61eea001cd1fb.tar.gz
    )
    FetchContent_GetProperties(stattools)
    if(NOT stattools_POPULATED)
        FetchContent_Populate(stattools)
        # Print debugging stuff
        message(STATUS "stattools_SOURCE_DIR = ${stattools_SOURCE_DIR}")
        message(STATUS "stattools_BINARY_DIR = ${stattools_BINARY_DIR}")
        # Copy into libs/ so it gets vendored for CRAN later
        file(COPY ${stattools_SOURCE_DIR}/ DESTINATION ${STATTOOLS_DIR})
    endif()
    add_subdirectory(${STATTOOLS_DIR})
endif()

get_property(_allTargets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "All targets at this point (after stattools): ${_allTargets}")

# Fail now if statools was not found!
if(NOT TARGET stattools)
    message(FATAL_ERROR "stattools target was not created! Check FetchContent / vendoring setup.")
endif()

#######################
# Create core library
#######################

# assemble list of cpp and h files in core folder
file(GLOB_RECURSE CORE_SOURCES core/*.h core/*.cpp)

#Core library
add_library(coreLib STATIC ${CORE_SOURCES})

target_sources(coreLib PRIVATE ${CORE_SOURCES})
target_include_directories(coreLib PRIVATE core SYSTEM INTERFACE core)
target_link_libraries(coreLib PUBLIC coretools stattools)
target_compile_features(coreLib PUBLIC cxx_std_17)
set_target_properties(coreLib PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(coreLib PRIVATE -Wall -Wextra)

#######################
# Executable
#######################

# add executable to be built
add_executable(${PROJECT_NAME}
        main.cpp
        )

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
        coreLib
        )
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

