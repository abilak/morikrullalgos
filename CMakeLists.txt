cmake_minimum_required(VERSION 3.20)
project(morikrull LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(morikrull_lib src/morikrull.cpp)
target_include_directories(morikrull_lib PUBLIC include)

set(NORMALIZ_ROOT "" CACHE PATH "Normaliz installation prefix")
find_path(NORMALIZ_INCLUDE_DIR libnormaliz/cone.h HINTS ${NORMALIZ_ROOT}/include REQUIRED)
find_library(NORMALIZ_LIBRARY normaliz HINTS ${NORMALIZ_ROOT}/lib REQUIRED)
get_filename_component(NORMALIZ_LIBRARY_DIR ${NORMALIZ_LIBRARY} DIRECTORY)
find_library(FLINT_LIBRARY flint HINTS ${NORMALIZ_LIBRARY_DIR} REQUIRED)
find_library(GMPXX_LIBRARY gmpxx HINTS ${NORMALIZ_LIBRARY_DIR} REQUIRED)
find_library(GMP_LIBRARY gmp HINTS ${NORMALIZ_LIBRARY_DIR} REQUIRED)
find_library(EANTICXX_LIBRARY eanticxx HINTS ${NORMALIZ_LIBRARY_DIR})
find_library(EANTIC_LIBRARY eantic HINTS ${NORMALIZ_LIBRARY_DIR})

target_include_directories(morikrull_lib SYSTEM PRIVATE ${NORMALIZ_INCLUDE_DIR})
set(MATH_LIBRARIES ${NORMALIZ_LIBRARY})
if(EANTICXX_LIBRARY AND EANTIC_LIBRARY)
  list(APPEND MATH_LIBRARIES ${EANTICXX_LIBRARY} ${EANTIC_LIBRARY})
endif()
list(APPEND MATH_LIBRARIES ${FLINT_LIBRARY} ${GMPXX_LIBRARY} ${GMP_LIBRARY})
target_link_libraries(morikrull_lib PRIVATE ${MATH_LIBRARIES})

add_executable(morikrull src/main.cpp)
target_link_libraries(morikrull PRIVATE morikrull_lib)
set_target_properties(morikrull PROPERTIES BUILD_RPATH ${NORMALIZ_LIBRARY_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(morikrull_lib PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_test(NAME regression COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression.sh $<TARGET_FILE:morikrull>)
