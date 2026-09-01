#pragma once

#include <string>
#include <TTree.h>

namespace brill {

struct DssdMatchEvent {
	int num = 0;
	int flag[8] = {0};
	int front_strip[8] = {0};
	int back_strip[8] = {0};
	double energy[8] = {0.0};
	double time[8] = {0.0};
	double x[8] = {0.0};
	double y[8] = {0.0};
	double z[8] = {0.0};
	// 0-f1b1, 1-f1b2, 2-f2b1, 3-f2b2
	int merge_tag[8] = {0};
	int run;
	int entry;
};

void SetupInput(TTree *tree, DssdMatchEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, DssdMatchEvent &event);
void Reset(DssdMatchEvent &event);

} // namespace brill