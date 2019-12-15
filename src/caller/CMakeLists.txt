file(GLOB SOURCES "*.cc")

add_library(caller
        STATIC
#        SHARED
        ${SOURCES}
)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory("../vcf" ${CMAKE_CURRENT_SOURCE_DIR}/vcf)
add_subdirectory("../fasta" ${CMAKE_CURRENT_SOURCE_DIR}/fasta)


if(APPLE)
    target_link_libraries(
        caller
        vcf
        fasta
        ${HTSLIB_LIBRARIES}
        ${ZLIB_LIBRARIES}
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
        ${BZIP2_LIBRARIES}
        ${TBB_IMPORTED_TARGETS}
    )
endif()
