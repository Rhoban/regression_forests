cmake_minimum_required(VERSION 3.16.3)
project(rhoban_regression_forests)
# Build options option(RHOBAN_RF_USES_GP "Is support for gaussian processes activated" OFF) option(RHOBAN_RF_BUILD_VIEWER "Build regression forests Viewer (require SFML)" OFF)
option(RHOBAN_RF_BUILD_EXAMPLES "Build examples/tools" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# Directories containing sources files
set(DIRECTORIES
  src/rhoban_regression_forests/algorithms
  src/rhoban_regression_forests/approximations
  src/rhoban_regression_forests/core
  src/rhoban_regression_forests/tools
)

# Build ALL_SOURCES
set (SOURCES)
set (ALL_SOURCES)
set (PREFIXED_SOURCES)
foreach (DIRECTORY ${DIRECTORIES})
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/${DIRECTORY}")
    include (${DIRECTORY}/sources.cmake)

    foreach (SOURCE ${SOURCES})
        set (PREFIXED_SOURCES ${PREFIXED_SOURCES} ${DIRECTORY}/${SOURCE})
    endforeach (SOURCE)

    set (ALL_SOURCES ${ALL_SOURCES} ${PREFIXED_SOURCES})
endforeach (DIRECTORY)

# if rhoban_viewer is found, then add the support for viewers
if (BUILD_VIEWER)
  set (ALL_SOURCES
    ${ALL_SOURCES}
    src/rhoban_regression_forests/tools/viewer.cpp
    )
endif(BUILD_VIEWER)

# Declare a C++ library
add_library(rhoban_regression_forests
  ${ALL_SOURCES}
)
target_link_libraries(rhoban_regression_forests
  rhoban_utils
  rhoban_random
)
target_include_directories(rhoban_regression_forests PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if (RHOBAN_RF_USES_GP)
  target_link_libraries(rhoban_regression_forests PUBLIC rhoban_gp)
  target_compile_definitions(rhoban_regression_forests PUBLIC RHOBAN_RF_USES_GP)
endif()

if(RHOBAN_RF_BUILD_VIEWER)
  target_link_libraries(rhoban_regression_forests PUBLIC rhoban_viewer)
endif()


# Testing input output functions
if(RHOBAN_RF_BUILD_EXAMPLES)
  add_executable(test_rf_io src/test_io.cpp)
  target_link_libraries(test_rf_io rhoban_regression_forests)
endif()

# if rhoban_viewer is found, then add the viewer executable as target
if (RHOBAN_RF_BUILD_VIEWER)
  add_executable(viewer src/viewer.cpp)
  target_link_libraries(viewer PRIVATE rhoban_regression_forests)
endif()
