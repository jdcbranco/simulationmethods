cmake_minimum_required(VERSION 3.3)
project(simulationmethods)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -${CMAKE_CXX_FLAGS}-std=c++11)
add_executable(main main.cpp Random.cpp Normal.cpp Option.h Vanilla.h Model.h Path.h MCModel.h Simulator.h ModelParams.h ModelResult.h BSVanilla.h Asian.h BSAsian.h SensitivityModel.h)
