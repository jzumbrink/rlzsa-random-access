cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)

# sdsl
set(CMAKE_WARN_DEPRECATED OFF CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/sdsl-lite/")
set_target_properties(
  uninstall-sdsl uninstall divsufsort divsufsort64
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1
)
target_include_directories(sdsl PUBLIC "${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/")

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++17")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

# ips4o
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ips4o")
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# unordered_dense
add_subdirectory("${CMAKE_SOURCE_DIR}/external/unordered_dense/")

# sais-lite-lcp
add_library(sais-lite-lcp INTERFACE)
target_include_directories(sais-lite-lcp INTERFACE "${CMAKE_SOURCE_DIR}/external/sais-lite-lcp/")

include_directories("${CMAKE_SOURCE_DIR}/external/ips4o/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/unordered_dense/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/sais-lite-lcp/")

add_executable(ri-build ri-build.cpp)
TARGET_LINK_LIBRARIES(ri-build sdsl)
TARGET_LINK_LIBRARIES(ri-build divsufsort)
TARGET_LINK_LIBRARIES(ri-build divsufsort64)
TARGET_LINK_LIBRARIES(ri-build ips4o)
TARGET_LINK_LIBRARIES(ri-build unordered_dense)
TARGET_LINK_LIBRARIES(ri-build sais-lite-lcp)

add_executable(ri-locate ri-locate.cpp)
TARGET_LINK_LIBRARIES(ri-locate sdsl)
TARGET_LINK_LIBRARIES(ri-locate divsufsort)
TARGET_LINK_LIBRARIES(ri-locate divsufsort64)

add_executable(ri-count ri-count.cpp)
TARGET_LINK_LIBRARIES(ri-count sdsl)
TARGET_LINK_LIBRARIES(ri-count divsufsort)
TARGET_LINK_LIBRARIES(ri-count divsufsort64)

#add_executable(ri-space ri-space.cpp)
#TARGET_LINK_LIBRARIES(ri-space sdsl)
#TARGET_LINK_LIBRARIES(ri-space divsufsort)
#TARGET_LINK_LIBRARIES(ri-space divsufsort64)

enable_testing()
include(GoogleTest)
include(Dart)

add_executable(test "${CMAKE_SOURCE_DIR}/tests/test.cpp")

TARGET_LINK_LIBRARIES(test sdsl)
TARGET_LINK_LIBRARIES(test divsufsort)
TARGET_LINK_LIBRARIES(test divsufsort64)
TARGET_LINK_LIBRARIES(test gtest_main)

set_target_properties(test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/")

gtest_discover_tests(test)