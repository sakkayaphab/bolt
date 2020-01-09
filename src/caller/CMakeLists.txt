file(GLOB SOURCES "*.cc")

add_library(caller
        STATIC
#        SHARED
        ${SOURCES}
)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
#message("PROJECT_SOURCE_DIR : " ${PROJECT_SOURCE_DIR}/src)
#add_subdirectory("../vcf" ${CMAKE_CURRENT_SOURCE_DIR}/vcf)
#add_subdirectory("../fasta" ${CMAKE_CURRENT_SOURCE_DIR}/fasta)

#find_library(HTSLIB_LIBRARIES hts HINTS ${LIBRARY_LINK_PREFIX})
#find_library(TBB_IMPORTED_TARGETS tbb HINTS ${LIBRARY_LINK_PREFIX})
#find_library(ZLIB_LIBRARIES z HINTS ${LIBRARY_LINK_PREFIX})
#find_library(BZIP2_LIBRARIES bz2 HINTS ${LIBRARY_LINK_PREFIX})

if(APPLE)
    target_link_libraries(
        caller
        vcf
        fasta
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        #${BZIP2_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}
    )
elseif(UNIX AND NOT APPLE)
    target_link_libraries(
        caller
        vcf
        fasta
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        #${BZIP2_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}
    )
endif()
