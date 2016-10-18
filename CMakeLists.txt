cmake_minimum_required(VERSION 3.6)
project(Coiote_heuristic)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp heuristic.h heuristic.cpp utils.h)
add_executable(Coiote_heuristic ${SOURCE_FILES})