cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 17)

#set(CMAKE_VERBOSE_MAKEFILE on)
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_MACOSX_RPATH 1)
#include(GNUInstallDirs)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

# add external project
#include(ExternalProject)

#set(CMAKE_CXX_STANDARD 14)
set(LIBRARY_LINK_PREFIX "/usr/local/lib/")
#set(LIBRARY_LINK_PREFIX2 "/usr/lib/")
find_library(HTSLIB_LIBRARIES hts HINTS ${LIBRARY_LINK_PREFIX})
find_library(TBB_IMPORTED_TARGETS tbb HINTS ${LIBRARY_LINK_PREFIX})
find_library(ZLIB_LIBRARIES z HINTS ${LIBRARY_LINK_PREFIX})
find_library(BZIP2_LIBRARIES bz2 HINTS ${LIBRARY_LINK_PREFIX})

#link_directories(${LIBRARY_LINK_PREFIX})
message("LIBRARY_LINK_PREFIX : " ${LIBRARY_LINK_PREFIX})
#set(INCLUDE_LIBRARY_PREFIX "/usr/local/include/")
#set(INCLUDE_LIBRARY_PREFIX2 "/usr/include/")
include_directories(${INCLUDE_LIBRARY_PREFIX})

message("INCLUDE_LIBRARY_PREFIX : " ${INCLUDE_LIBRARY_PREFIX})
message(STATUS $ENV{$LD_LIBRARY_PATH})

set(CMAKE_INSTALL_RPATH ${LIBRARY_LINK_PREFIX})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_CXX_FLAGS "-O3")
#add_executable(jinn main.cpp)

# Build external library
#include(cmake/FindHTSlib.cmake)
#include(cmake/FindZLIBlib.cmake)
#include(cmake/FindLZMAlib.cmake)
#include(cmake/FindBZip2.cmake)
#include(cmake/FindTBB.cmake)

#enable_testing()

#include header
#include_directories(${HTSLIB_INCLUDE_DIRS})

# build bolt
add_subdirectory(src)

#option(BUILD_TESTS "Build test" OFF)
#if(BUILD_TESTS)
#    add_subdirectory(test)
#endif(BUILD_TESTS)

