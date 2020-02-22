file(GLOB SOURCES "*.cc")

add_executable(bolt ${SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(caller)
add_subdirectory(fasta)
add_subdirectory(graph)

target_link_libraries(bolt
        caller
        graph
  )

if (NOT ${INSTALL_BIN_PREFIX} STREQUAL "")
    message("INSTALL BINARY PATH : " ${INSTALL_BIN_PREFIX})
    set(INSTALL_BIN_DIR ${INSTALL_BIN_PREFIX} CACHE PATH "Installation directory for executables")
    install(TARGETS bolt DESTINATION ${INSTALL_BIN_DIR})
else()
    message("INSTALL BINARY PATH : " "/usr/local/bin")
    set(INSTALL_BIN_DIR "/usr/local/bin" CACHE PATH "Installation directory for executables")
    install(TARGETS bolt DESTINATION ${INSTALL_BIN_DIR})
endif()