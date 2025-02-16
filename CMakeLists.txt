cmake_minimum_required(VERSION 2.6)
set(CMAKE_CXX_STANDARD 20)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (s-rlbwt)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/internal)

# malloc_count
add_library(malloc_count INTERFACE)
target_sources(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/malloc_count.c")
target_include_directories(malloc_count INTERFACE "${CMAKE_SOURCE_DIR}/external/malloc_count/")
target_link_libraries(malloc_count INTERFACE dl)

# sdsl
include_directories("~/include")
link_directories("~/lib")

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -ggdb -fstrict-aliasing -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing")

# ips4o
add_subdirectory("${CMAKE_SOURCE_DIR}/external/ips4o")
set_target_properties(ips4o_example PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

# unordered_dense
add_subdirectory("${CMAKE_SOURCE_DIR}/external/unordered_dense/")

# libsais
set(LIBSAIS_USE_OPENMP OFF CACHE BOOL "" FORCE)
add_subdirectory("${CMAKE_SOURCE_DIR}/external/libsais/")

include_directories("${CMAKE_SOURCE_DIR}/external/ips4o/include/")
include_directories("${CMAKE_SOURCE_DIR}/external/unordered_dense/include/")

add_executable(rlzsa-build rlzsa-build.cpp)
TARGET_LINK_LIBRARIES(rlzsa-build sdsl)
TARGET_LINK_LIBRARIES(rlzsa-build divsufsort)
TARGET_LINK_LIBRARIES(rlzsa-build divsufsort64)
TARGET_LINK_LIBRARIES(rlzsa-build ips4o)
TARGET_LINK_LIBRARIES(rlzsa-build unordered_dense)
TARGET_LINK_LIBRARIES(rlzsa-build libsais)
TARGET_LINK_LIBRARIES(rlzsa-build malloc_count)

add_executable(rlzsa-ra rlzsa-random-access.cpp)
TARGET_LINK_LIBRARIES(rlzsa-ra sdsl)
TARGET_LINK_LIBRARIES(rlzsa-ra divsufsort)
TARGET_LINK_LIBRARIES(rlzsa-ra divsufsort64)
TARGET_LINK_LIBRARIES(rlzsa-ra ips4o)
TARGET_LINK_LIBRARIES(rlzsa-ra unordered_dense)
TARGET_LINK_LIBRARIES(rlzsa-ra libsais)
TARGET_LINK_LIBRARIES(rlzsa-ra malloc_count)