cmake_minimum_required(VERSION "3.14")

project(lsoda LANGUAGES Fortran CXX)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

set(CMAKE_CXX_STANDARD 11)

# lsoda
add_library(lsoda SHARED 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LSODA.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/wrapper.cpp
)
set_target_properties(lsoda PROPERTIES PREFIX "lib")

# dop853
add_library(dop853 SHARED
  src/dop853_constants.f90
  src/dop853_module.f90
  src/dop853_c_interface.f90
)
set_target_properties(dop853 PROPERTIES PREFIX "lib")

if (SKBUILD)
  install(TARGETS lsoda DESTINATION numbalsoda)
  install(TARGETS dop853 DESTINATION numbalsoda)
else()
  install(TARGETS lsoda DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/numbalsoda/)
  install(TARGETS dop853 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/numbalsoda/)
endif()