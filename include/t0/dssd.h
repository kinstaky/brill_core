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
	DssdNormalizeParameters(const int fs, const int bs);
	inline double NormEnergy(
		const int side, const int strip, const double raw_energy
	) const {
		return side == 0
			? front_p0[strip]
				+ front_p1[strip] * raw_energy
				+ front_p2[strip] * raw_energy * raw_energy
				+ front_p3[strip] * raw_energy * raw_energy * raw_energy
			: back_p0[strip]
				+ back_p1[strip] * raw_energy
				+ back_p2[strip] * raw_energy * raw_energy
				+ back_p3[strip] * raw_energy * raw_energy * raw_energy;
	}
	int Write(const std::string &path) const;
	int Read(const std::string &path);
	void Apply(const DssdEvent &input, DssdEvent &output) const;

	int front_strips = 0;
	int back_strips = 0;
	double front_p0[kMaxStrips] = {0.0};
	double front_p1[kMaxStrips] = {1.0};
	double front_p2[kMaxStrips] = {0.0};
	double front_p3[kMaxStrips] = {0.0};
	double back_p0[kMaxStrips] = {0.0};
	double back_p1[kMaxStrips] = {1.0};
	double back_p2[kMaxStrips] = {0.0};
	double back_p3[kMaxStrips] = {0.0};
};

class CalibrationParameters {
public:
	CalibrationParameters(int layers);

	inline size_t Layers() const { return p0.size(); }
	int Write(const std::string &path) const;
	int Read(const std::string &path);

	std::vector<double> p0;
	std::vector<double> p1;
};


void MatchDssdEvent(
	const DssdEvent &input,
	const SiliconDetectorConfig &detector,
	DssdMatchEvent &output
);

} // namespace brill
