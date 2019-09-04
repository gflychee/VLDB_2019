## Main build

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/git_info.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp" @ONLY)

set( HEADER_FILES
  git_info.hpp
  prelude.hpp
  termcolor.hpp
  logging.hpp
  rand.hpp
  experiment_reporter.hpp
  io.hpp
  connected_components.hpp
  km_sampler.hpp
  cluster_vertex.hpp
  guesser.hpp
  kmpc.hpp
  Alo4.hpp
  counts_cache.hpp
)

set( LIB_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp"
  io.cpp
  logging.cpp
  rand.cpp
  connected_components.cpp
  km_sampler.cpp
  scores.cpp
)
add_library(ugraph STATIC ${HEADER_FILES} ${LIB_FILES} )


add_executable(ugraph-Alo4 Alo4.cpp)
target_link_libraries(ugraph-Alo4
        ugraph
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_DATE_TIME_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        ${Boost_IOSTREAMS_LIBRARY}
        ${BZIP2_LIBRARIES})



