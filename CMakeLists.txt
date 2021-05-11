cmake_minimum_required(VERSION 3.19)
project(Sim)

message(STATUS "GCC Version - ${CMAKE_CXX_COMPILER_VERSION}")

set(CMAKE_CXX_STANDARD 17)
# BOOST
set(Boost_USE_MULTITHREADED ON)
set(BOOST_ROOT "/opt/boost/boost_1_76_0")
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
    message(FATAL_ERROR "Fatal error: Boost (version >= 1.55) required.")
else()
    message(STATUS "Setting up BOOST")
    message(STATUS " Includes - ${Boost_INCLUDE_DIRS}")
    message(STATUS " Library  - ${Boost_LIBRARY_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif (NOT Boost_FOUND)

# FMT
find_package(fmt)
if (NOT fmt_FOUND)
    message(FATAL_ERROR "Fatal error: fmt required.")
else()
    message(STATUS "Setting up fmt")
    include_directories(${fmt_INCLUDE_DIRS})
endif(NOT fmt_FOUND)

add_executable(Sim SimpleMD.cpp Engine.cpp Engine.h Particle.cpp Particle.h ParticleState.cpp ParticleState.h Vector.h Vector.cpp setup.h HexagonalGrid.h)