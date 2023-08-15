# This file is part of **mygcbgcmc**.
# 
# mygcbgcmc is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# mygcbgcmc is distributed with the hope that it will be useful,
# but **WITHOUT ANY WARRANTY**; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
# 
# See the GNU Lesser General Public License for more details.


cmake_minimum_required(VERSION 3.17)

project(cppgcmc CXX)

set(SOURCES charmm.hpp const.hpp gcmccall.hpp gcmcmole.hpp gcmcmove.hpp
    grid.hpp mygcbgcmc.hpp myrand.hpp charmm.cpp gcmccall.cpp gcmcmain.cpp
    gcmcmole.cpp gcmcmove.cpp grid.cpp mygcbgcmc.cpp MersenneTwister.h)


find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# make sure to build with -O3
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  add_compile_options(-O3)
endif()

add_executable(gcbgcmc.exe ${SOURCES})


