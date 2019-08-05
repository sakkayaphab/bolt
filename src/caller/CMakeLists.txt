file(GLOB SOURCES "*.cc")

include_directories(${CMAKE_BINARY_DIR}/src/)


add_library(caller
#        STATIC
        SHARED
    ${SOURCES}
)

add_subdirectory("../vcf" ${CMAKE_CURRENT_BINARY_DIR}/vcf)
add_subdirectory("../fasta" ${CMAKE_CURRENT_BINARY_DIR}/fasta)

if(APPLE)
    target_link_libraries(
        caller
        vcf
        fasta
        
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${BZIP2_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}

    )
elseif(UNIX AND NOT APPLE)
    target_link_libraries(
        caller
        vcf
        fasta
        
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
        ${BZIP2_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}

    )
endif()
