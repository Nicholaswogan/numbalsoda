cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Release)

add_definitions( -std=c++11 ) 
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(lsoda SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/LSODA.cpp
                         ${CMAKE_CURRENT_SOURCE_DIR}/src/wrapper.cpp)

set_target_properties(lsoda PROPERTIES SUFFIX ".so")

install(TARGETS lsoda DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/NumbaLSODA/)




