cmake_minimum_required(VERSION 3.17)
project(violin_boltzman)

set(CMAKE_CXX_STANDARD 11)

add_executable(violin_boltzman LatticeBoltzmann2D.cpp LatticeBoltzmann.cpp LatticeBoltzmann3D.cpp)