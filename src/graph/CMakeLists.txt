file(GLOB SOURCES "*.cc")

add_library(graph
        STATIC
        ${SOURCES}
        )