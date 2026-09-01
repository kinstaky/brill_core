#ifndef __T0_EVENT_H__
#define __T0_EVENT_H__

#include <string>
#include <TTree.h>

namespace brill {

struct T0Event {
	int num;
	int layer[8];
	int flag[8];
	int charge[8];
	int mass[8];
	double energy[8][6];
	double time[8][6];
	double x[8][4];
	double y[8][4];
	double z[8][4];
	int last[8][4];
	int run;
	int entry;
};

void SetupInput(TTree *tree, T0Event &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, T0Event &event);
void Reset(T0Event &event);

}

#endif