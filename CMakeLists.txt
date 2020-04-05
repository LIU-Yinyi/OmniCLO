cmake_minimum_required(VERSION 3.10)
project(OmniCLO)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_BUILD_TYPE Release)

#add_definitions("-O3")

set(INC_DIR_DISPLAY display)
set(INC_DIR_PHYSICS physics)
set(INC_DIR_PROBLEM problems)
set(LIB_LIST pagmo pthread)

include_directories(/usr/local/include ./ ${INC_DIR_DISPLAY} ${INC_DIR_PHYSICS} ${INC_DIR_PROBLEM})
link_directories(/usr/local/lib)

add_executable(OmniCLO main.cpp)
target_link_libraries(OmniCLO ${LIB_LIST})

file(GLOB SRC_UNIT_TEST "test/*.cpp")
add_executable(UnitTest ${SRC_UNIT_TEST})
target_link_libraries(UnitTest ${LIB_LIST})
