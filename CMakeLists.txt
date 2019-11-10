cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 17)

#set(CMAKE_VERBOSE_MAKEFILE on)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH 1)
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# add external project
include(ExternalProject)

# Build external library
include(cmake/FindHTSlib.cmake)
include(cmake/FindZLIBlib.cmake)
include(cmake/FindLZMAlib.cmake)
include(cmake/FindBZip2.cmake)
include(cmake/FindTBB.cmake)

enable_testing()

#include header
include_directories(${HTSLIB_INCLUDE_DIRS})

# build bolt
add_subdirectory(src)

option(BUILD_TESTS "Build test" OFF)
if(BUILD_TESTS)
    add_subdirectory(test)
endif(BUILD_TESTS)
