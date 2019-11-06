project (bolt)

file(GLOB SOURCES "*.cc")

add_executable(bolt ${SOURCES})

include_directories(caller)

add_subdirectory(caller)
add_subdirectory(bwt)

target_link_libraries(bolt
  PUBLIC
    caller
  )