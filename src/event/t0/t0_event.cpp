#include "include/event/t0/t0_event.h"

namespace brill {

void SetupInput(TTree *tree, T0Event &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "layer").c_str(), event.layer);
	tree->SetBranchAddress((prefix + "flag").c_str(), event.flag);
	tree->SetBranchAddress((prefix + "charge").c_str(), event.charge);
	tree->SetBranchAddress((prefix + "mass").c_str(), event.mass);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "time").c_str(), event.time);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
	tree->SetBranchAddress((prefix + "z").c_str(), event.z);
	tree->SetBranchAddress((prefix + "last").c_str(), event.last);
	tree->SetBranchAddress((prefix + "run").c_str(), &event.run);
	tree->SetBranchAddress((prefix + "entry").c_str(), &event.entry);
}

void SetupOutput(TTree *tree, T0Event &event) {
	tree->Branch("num", &event.num, "num/I");
	tree->Branch("layer", event.layer, "l[num]/I");
	tree->Branch("flag", event.flag, "flag[num]/I");
	tree->Branch("charge", event.charge, "Z[num]/I");
	tree->Branch("mass", event.mass, "A[num]/I");
	tree->Branch("energy", event.energy, "e[num][6]/D");
	tree->Branch("time", event.time, "t[num][6]/D");
	tree->Branch("x", event.x, "x[num][4]/D");
	tree->Branch("y", event.y, "y[num][4]/D");
	tree->Branch("z", event.z, "z[num][4]/D");
	tree->Branch("last", event.last, "last[num][4]/I");
	tree->Branch("run", &event.run, "run/I");
	tree->Branch("entry", &event.entry, "entry/I");
}

void Reset(T0Event &event) {
	event.num = 0;
	for (int i = 0; i < 8; ++i) {
		event.layer[i] = 0;
		event.flag[i] = 0;
		event.charge[i] = 0;
		event.mass[i] = 0;
		for (int j = 0; j < 6; ++j) {
			event.energy[i][j] = 0.0;
			event.time[i][j] = 0.0;
		}
		for (int j = 0; j < 4; ++j) {
			event.x[i][j] = 0.0;
			event.y[i][j] = 0.0;
			event.z[i][j] = 0.0;
			event.last[i][j] = -1;
		}
	}
	event.run = -1;
	event.entry = -1;
}

} // namespace brill
