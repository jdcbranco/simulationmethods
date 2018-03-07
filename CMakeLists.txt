cmake_minimum_required(VERSION 3.9)
project(simulationmethods)

set(CMAKE_CXX_STANDARD 14)

add_executable(main main.cpp Random.cpp Normal.cpp Option.h Vanilla.h Model.cpp Model.h Path.h MCModel.cpp MCModel.h Simulator.h ModelParams.h ModelResult.h BSVanilla.h Asian.h)
add_executable(european european.cpp Random.cpp Normal.cpp Option.h Vanilla.h Model.cpp Model.h Path.h MCModel.cpp MCModel.h Simulator.h ModelParams.h ModelResult.h BSVanilla.h Asian.h)
add_executable(asian asian.cpp Random.cpp Normal.cpp Option.h Vanilla.h Model.cpp Model.h Path.h MCModel.cpp MCModel.h Simulator.h ModelParams.h ModelResult.h BSVanilla.h Asian.h)
add_executable(american asian.cpp Random.cpp Normal.cpp Option.h Vanilla.h Model.cpp Model.h Path.h MCModel.cpp MCModel.h Simulator.h ModelParams.h ModelResult.h BSVanilla.h Asian.h)