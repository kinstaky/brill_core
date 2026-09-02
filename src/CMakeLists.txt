add_subdirectory(event)
add_subdirectory(energy_calculator)
add_subdirectory(rebuild)
add_subdirectory(t0)

add_library(brill_config config.cpp)
target_include_directories(brill_config PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})

add_library(brill_core INTERFACE)
target_link_libraries(brill_core INTERFACE brill_config)

add_library(brill_particle INTERFACE)
target_link_libraries(brill_particle INTERFACE brill_particle_event)