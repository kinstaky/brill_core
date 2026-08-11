#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/t0/dssd_match_event.h"

namespace brill {

class DssdNormalizeParameters {
public:
	DssdNormalizeParameters(const int front_strip, const int back_strip);
	inline double NormEnergy(
		const int side, const int strip, const double raw_energy
	) const {
		return side == 0
			? front_p0_[strip] + front_p1_[strip] * raw_energy
			: back_p0_[strip] + back_p1_[strip] * raw_energy;
	}
	int Write(const std::string &path) const;
	int Read(const std::string &path);
	void Apply(const DssdEvent &input, DssdEvent &output) const;

	int front_strips_ = 0;
	int back_strips_ = 0;
	double front_p0_[kMaxStrips] = {0.0};
	double front_p1_[kMaxStrips] = {1.0};
	double front_p2_[kMaxStrips] = {0.0};
	double back_p0_[kMaxStrips] = {0.0};
	double back_p1_[kMaxStrips] = {1.0};
	double back_p2_[kMaxStrips] = {0.0};
};


void MatchDssdEvent(
	const DssdEvent &input,
	const SiliconDetectorConfig &detector,
	DssdMatchEvent &output
);

} // namespace brill
