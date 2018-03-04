cmake_minimum_required(VERSION 3.9)
project(simulationmethods)

set(CMAKE_CXX_STANDARD 14)

add_executable(european european.cpp)
add_executable(asian asian.cpp)