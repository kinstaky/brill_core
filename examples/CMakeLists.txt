add_executable(match_dssd match_dssd.cpp)
target_link_libraries(
	match_dssd
	PRIVATE
	brill_config
	brill_t0_dssd
)

add_executable(track_t0 track_t0.cpp)
target_link_libraries(
	track_t0
	PRIVATE
	brill_config
	brill_dssd_match_event
	brill_t0_event
	brill_silicon_event
	ROOT::Graf
)

add_executable(rebuild_t0 rebuild_t0.cpp)
target_link_libraries(
	rebuild_t0
	PRIVATE
	brill_config
	brill_t0_event
	brill_particle_event
	brill_energy_calculator
)

add_executable(calibrate_t0 calibrate_t0.cpp)
target_link_libraries(
	calibrate_t0
	PRIVATE
	brill_config
	brill_t0_event
	brill_energy_calculator
)

add_subdirectory(estimate)
# add_subdirectory(rebuild)