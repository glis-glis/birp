cmake_minimum_required(VERSION 3.14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Policies (for properly setting the timestamps of extracted contents)
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)
if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif (POLICY CMP0135)

# To avoid warnings when compiling on MacOS on CRAN
set(CMAKE_OSX_DEPLOYMENT_TARGET "13.7" CACHE STRING "Minimum OS X deployment version")

# Project
project(birp LANGUAGES CXX)

# Download coretools and stattools
include(FetchContent)
FetchContent_Declare(coretools
        URL https://bitbucket.org/wegmannlab/coretools/get/bdb81499c14964952341da8fdefe547cb69634b6.tar.gz
        SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/coretools")

FetchContent_Declare(stattools
        URL https://bitbucket.org/wegmannlab/stattools/get/2e11bbc29def51fac86d178e4cdaa4e65248615d.tar.gz
        SOURCE_DIR  "${CMAKE_CURRENT_SOURCE_DIR}/stattools")

FetchContent_MakeAvailable(coretools stattools)

#######################
# Create core library
#######################

# assemble list of cpp and h files in core folder
file(GLOB_RECURSE CORE_SOURCES core/*.h core/*.cpp)

#Core library
add_library(${PROJECT_NAME}Lib STATIC ${CORE_SOURCES})

target_sources(${PROJECT_NAME}Lib PRIVATE ${CORE_SOURCES})
target_include_directories(${PROJECT_NAME}Lib PRIVATE core SYSTEM INTERFACE core)
target_link_libraries(${PROJECT_NAME}Lib PUBLIC coretools stattools)
target_compile_features(${PROJECT_NAME}Lib PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME}Lib PROPERTIES CXX_EXTENSIONS OFF)
target_compile_options(${PROJECT_NAME}Lib PRIVATE -Wall -Wextra)

#######################
# Executable
#######################

# add executable to be built
add_executable(${PROJECT_NAME} main.cpp)

# link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${PROJECT_NAME}Lib)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)

