cmake_minimum_required(VERSION 3.4.0)

project(lsoda)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(CMAKE_CXX_STANDARD 11)
add_library(lsoda SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/LSODA.cpp
                         ${CMAKE_CURRENT_SOURCE_DIR}/src/wrapper.cpp)
set_target_properties(lsoda PROPERTIES PREFIX "lib")

if (SKBUILD)
  install(TARGETS lsoda DESTINATION NumbaLSODA)
else()
  install(TARGETS lsoda DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/NumbaLSODA/)
endif()