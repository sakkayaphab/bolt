project (bolt)

file(GLOB SOURCES "*.cc")

add_executable(bolt ${SOURCES})

include_directories(caller)

add_subdirectory(caller)
#add_subdirectory(bwt)

target_link_libraries(bolt
  PUBLIC
    caller
  )

if (NOT ${INSTALL_PREFIX} STREQUAL "")
    message("INSTALL BINARY PATH : " ${INSTALL_PREFIX})
    set(INSTALL_BIN_DIR ${INSTALL_BIN_PREFIX} CACHE PATH "Installation directory for executables")
    install(TARGETS bolt DESTINATION ${INSTALL_BIN_DIR})
else()
    message("INSTALL BINARY PATH : " INSTALL_PREFIX)
    set(INSTALL_BIN_DIR INSTALL_PREFIX CACHE PATH "Installation directory for executables")
    install(TARGETS bolt DESTINATION ${INSTALL_BIN_DIR})
endif()