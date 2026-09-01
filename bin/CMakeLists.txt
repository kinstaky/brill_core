add_executable(normalize normalize.cpp)
target_link_libraries(
	normalize
	PRIVATE
	brill_config
	brill_dssd_event
	brill_t0_dssd
	ROOT::Graf
)

add_subdirectory(estimate)