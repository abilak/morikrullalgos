cmake_minimum_required(VERSION 3.20)
project(morikrull LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(morikrull_lib src/morikrull.cpp)
target_include_directories(morikrull_lib PUBLIC include)

find_path(BOOST_INCLUDE_DIR boost/multiprecision/cpp_int.hpp
  PATHS /opt/homebrew/include /usr/local/include)
if(NOT BOOST_INCLUDE_DIR)
  message(FATAL_ERROR "Boost.Multiprecision headers are required")
endif()
target_include_directories(morikrull_lib SYSTEM PRIVATE ${BOOST_INCLUDE_DIR})

add_executable(morikrull src/main.cpp)
target_link_libraries(morikrull PRIVATE morikrull_lib)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(morikrull_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME regression COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.sh $<TARGET_FILE:morikrull>)
