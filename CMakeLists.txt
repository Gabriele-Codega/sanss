cmake_minimum_required(VERSION 3.13.4)

project(sanss)

find_package(deal.II 9.6.0 REQUIRED
    HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
)
deal_ii_initialize_cached_variables()

include_directories(include)

set(SOURCES
     main.cc
     source/navier_stokes.cc
     source/time_control.cc
    )


add_executable(sanss ${SOURCES})
deal_ii_setup_target(sanss)

set_target_properties(sanss PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
