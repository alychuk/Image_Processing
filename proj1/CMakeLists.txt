cmake_minimum_required(VERSION 3.13)  # CMake version check
project(proj1)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 14)            # Enable c++14 standard

add_library( STATIC)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(proj1 driver.cpp)
