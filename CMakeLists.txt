cmake_minimum_required(VERSION 3.7)
project(Linux)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        ANNOY/annoylib.h
        ANNOY/kissrandom.h
        LargeVis.cpp
        LargeVis.h
        main.cpp)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -pthread -lgsl -lgslcblas -Ofast -march=native -ffast-math")
include_directories(/usr/local/include)
link_directories(/usr/local/lib)
add_executable(Linux ${SOURCE_FILES})