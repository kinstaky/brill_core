add_executable(estimate_normalize estimate_normalize.cpp)
target_link_libraries(
	estimate_normalize
	PRIVATE
	brill_config
	brill_t0_dssd
	ROOT::Graf
)
