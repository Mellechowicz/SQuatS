cmake_minimum_required(VERSION 3.20)
project(exsqs VERSION 0.0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(exsqs src/main.cpp)
target_include_directories(exsqs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
