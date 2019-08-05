file(GLOB SOURCES "*.cc")

add_library(bwt
#        STATIC
        SHARED
    ${SOURCES}
)