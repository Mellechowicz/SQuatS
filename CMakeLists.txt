cmake_minimum_required(VERSION 3.20)
project(exsqs VERSION 0.2.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)

set(SPGLIB_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_Fortran OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_Python OFF CACHE BOOL "" FORCE)
set(SPGLIB_WITH_TESTS OFF CACHE BOOL "" FORCE)
set(SPGLIB_USE_OMP OFF CACHE BOOL "" FORCE)
FetchContent_Declare(spglib
  GIT_REPOSITORY https://github.com/spglib/spglib.git
  GIT_TAG v2.7.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(spglib)

add_library(exsqs_core
  src/structure.cpp
  src/lattice.cpp
  src/zones.cpp
  src/symmetry.cpp
  src/displacements.cpp
  src/dedup.cpp
)
target_include_directories(exsqs_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(TARGET Spglib::symspg)
  target_link_libraries(exsqs_core PUBLIC Spglib::symspg)
else()
  target_link_libraries(exsqs_core PUBLIC symspg)
endif()

add_executable(exsqs src/main.cpp)
target_link_libraries(exsqs PRIVATE exsqs_core)
