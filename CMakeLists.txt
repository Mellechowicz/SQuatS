cmake_minimum_required(VERSION 3.20)
project(exsqs VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_library(exsqs_core
  src/structure.cpp
  src/lattice.cpp
)
target_include_directories(exsqs_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(exsqs src/main.cpp)
target_link_libraries(exsqs PRIVATE exsqs_core)
