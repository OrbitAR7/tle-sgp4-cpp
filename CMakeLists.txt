cmake_minimum_required(VERSION 3.14)
project(tle_sgp4_cpp 
    VERSION 1.0.0
    DESCRIPTION "C++ SGP4/SDP4 Orbit Propagator for TLE Data"
    LANGUAGES CXX)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
    add_compile_options(/W4)
endif()

# Header-only library interface
add_library(sgp4_propagator INTERFACE)
target_include_directories(sgp4_propagator INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Examples
option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    add_executable(iss_tracking examples/iss_tracking.cpp)
    target_link_libraries(iss_tracking PRIVATE sgp4_propagator)

    add_executable(constellation_demo examples/constellation_demo.cpp)
    target_link_libraries(constellation_demo PRIVATE sgp4_propagator)
endif()

# Tests
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_executable(test_sgp4 tests/test_sgp4.cpp)
    target_link_libraries(test_sgp4 PRIVATE sgp4_propagator)
    add_test(NAME SGP4Tests COMMAND test_sgp4)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS sgp4_propagator
    EXPORT sgp4_propagator_targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES include/sgp4_propagator.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
