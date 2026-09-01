add_executable(estimate_t0_center estimate_t0_center.cpp)
target_link_libraries(
	estimate_t0_center
	PRIVATE
	brill_config
	brill_t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_pid estimate_t0_pid.cpp)
target_link_libraries(
	estimate_t0_pid
	PRIVATE
	brill_config
	brill_silicon_event
	brill_t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_straight estimate_t0_straight.cpp)
target_link_libraries(
	estimate_t0_straight
	PRIVATE
	brill_config
	brill_silicon_event
	brill_t0_dssd
	ROOT::Graf
)

add_executable(estimate_t0_csi_pid estimate_t0_csi_pid.cpp)
target_link_libraries(
	estimate_t0_csi_pid
	PRIVATE
	brill_config
	brill_silicon_event
	ROOT::Graf
)
