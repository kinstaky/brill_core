#include "include/event/t0/dssd_match_event.h"

namespace brill {

void SetupInput(TTree *tree, DssdMatchEvent &event, const std::string &prefix) {
	tree->SetBranchAddress((prefix + "num").c_str(), &event.num);
	tree->SetBranchAddress((prefix + "flag").c_str(), event.flag);
	tree->SetBranchAddress((prefix + "front_strip").c_str(), event.front_strip);
	tree->SetBranchAddress((prefix + "back_strip").c_str(), event.back_strip);
	tree->SetBranchAddress((prefix + "energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix + "time").c_str(), event.time);
	tree->SetBranchAddress((prefix + "x").c_str(), event.x);
	tree->SetBranchAddress((prefix + "y").c_str(), event.y);
	tree->SetBranchAddress((prefix + "z").c_str(), event.z);
	tree->SetBranchAddress((prefix + "merge_tag").c_str(), event.merge_tag);
}

void SetupOutput(TTree *tree, DssdMatchEvent &event) {
	tree->Branch("num", &event.num, "num/I");
	tree->Branch("flag", event.flag, "flag[num]/I");
	tree->Branch("front_strip", event.front_strip, "fs[num]/I");
	tree->Branch("back_strip", event.back_strip, "bs[num]/I");
	tree->Branch("energy", event.energy, "e[num]/D");
	tree->Branch("time", event.time, "t[num]/D");
	tree->Branch("x", event.x, "x[num]/D");
	tree->Branch("y", event.y, "y[num]/D");
	tree->Branch("z", event.z, "z[num]/D");
	tree->Branch("merge_tag", event.merge_tag, "mt[num]/I");
}

void Reset(DssdMatchEvent &event) {
	event.num = 0;
}

} // namespace brill
