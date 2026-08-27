#pragma once

#include <string>

#include <TTree.h>

namespace brill {

struct DssdEvent {
	int front_num;
	int front_strip[8];
	double front_energy[8];
	double front_integral[8];
	double front_time[8];
	int back_num;
	int back_strip[8];
	double back_energy[8];
	double back_integral[8];
	double back_time[8];
	int run;
	int entry;
};

void SetupInput(TTree *tree, DssdEvent &event, const std::string &prefix = "");

void SetupOutput(TTree *tree, DssdEvent &event);

void Reset(DssdEvent &event);

void Update(DssdEvent &event, bool front, int strip, double energy, double time);

} // namespace glimmer
