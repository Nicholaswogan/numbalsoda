cmake_minimum_required(VERSION "3.14")

project(lsoda LANGUAGES Fortran CXX)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
set(CMAKE_CXX_STANDARD 11)

# lsoda
add_library(lsoda SHARED 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LSODA.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/wrapper.cpp
)
set_target_properties(lsoda PROPERTIES PREFIX "lib")

# dop853, which is used two times
add_library(dop853_mod
  src/dop853_constants.f90
  src/dop853_module.f90
)

# dop853
add_library(dop853 SHARED
  src/dop853_c_interface.f90
)
set_target_properties(dop853 PROPERTIES PREFIX "lib")
target_link_libraries(dop853 dop853_mod)

# solve_ivp
add_library(solve_ivp SHARED
  src/brent_mod.f90
  src/dop853_solve_ivp.f90
)
set_target_properties(solve_ivp PROPERTIES PREFIX "lib")
target_link_libraries(solve_ivp dop853_mod)

if (SKBUILD)
  install(TARGETS lsoda DESTINATION numbalsoda)
  install(TARGETS dop853 DESTINATION numbalsoda)
  install(TARGETS solve_ivp DESTINATION numbalsoda)
else()
  install(TARGETS lsoda DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/numbalsoda/)
  install(TARGETS dop853 DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/numbalsoda/)
  install(TARGETS solve_ivp DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/numbalsoda/)
endif()
