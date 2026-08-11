#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

bool InTrackWindow(
	const brill::DssdMatchEvent &left,
	int left_index,
	const brill::DssdMatchEvent &right,
	int right_index,
	const brill::TrackWindowConfig &window
) {
	double dx = right.x[right_index] - left.x[left_index];
	double dy = right.y[right_index] - left.y[left_index];
	return
		dx >= window.min
		&& dx <= window.max
		&& dy >= window.min
		&& dy <= window.max;
}

void FillPairPid(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	const brill::TrackWindowConfig &window,
	TH2F &histogram
) {
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			if (!InTrackWindow(left, i, right, j, window)) continue;
			histogram.Fill(right.energy[j], left.energy[i]);
		}
	}
}

void FillSiliconPid(
	const brill::DssdMatchEvent &dssd,
	const brill::SiliconEvent &silicon,
	TH2F &histogram
) {
	if (!silicon.valid) return;
	if (dssd.num == 1) {
		histogram.Fill(silicon.energy, dssd.energy[0]);
	}
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_pid", "Estimate T0 PID after DSSD matching.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		return 1;
	}
	if (result.count("trigger")) {
		config.root_.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const std::string match_dir = brill::JoinPath(config.root_.workspace, config.paths_.match);
	const std::string ingot_dir = brill::JoinPath(config.root_.workspace, config.paths_.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.root_.trigger);

	TChain chain1("tree");
	TChain chain2("tree");
	TChain chain3("tree");
	TChain chain4("tree");
	TChain chain_s("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
		++added_runs;
		chain1.Add(TString::Format(
			"%s/t0d1_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain2.Add(TString::Format(
			"%s/t0d2_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain3.Add(TString::Format(
			"%s/t0d3_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain4.Add(TString::Format(
			"%s/t0d4_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain_s.Add(TString::Format(
			"%s/t0s_%s%04d.root",
			ingot_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain3, "d3");
	chain1.AddFriend(&chain4, "d4");
	chain1.AddFriend(&chain_s, "s");

	brill::DssdMatchEvent event1;
	brill::DssdMatchEvent event2;
	brill::DssdMatchEvent event3;
	brill::DssdMatchEvent event4;
	brill::SiliconEvent event_s;
	brill::SetupInput(&chain1, event1);
	brill::SetupInput(&chain1, event2, "d2.");
	brill::SetupInput(&chain1, event3, "d3.");
	brill::SetupInput(&chain1, event4, "d4.");
	brill::SetupInput(&chain1, event_s, "s.");

	TString output_path = TString::Format(
		"%s/t0_pid_%s%04d_%04d.root",
		brill::JoinPath(config.root_.workspace, config.paths_.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TH2F d1d2_pid("d1d2", "D1-D2 PID", 5000, 0.0, 60000.0, 5000, 0.0, 30000.0);
	TH2F d2d3_pid("d2d3", "D2-D3 PID", 5000, 0.0, 50000.0, 5000, 0.0, 80000.0);
	TH2F d3d4_pid("d3d4", "D3-D4 PID", 5000, 0.0, 45000.0, 5000, 0.0, 45000.0);
	TH2F d4s_pid("d4s", "D4-S PID", 5000, 0.0, 80000.0, 5000, 0.0, 45000.0);

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Filling T0 PID   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		FillPairPid(event1, event2, config.track_.d2d1_window, d1d2_pid);
		FillPairPid(event2, event3, config.track_.d3d2_window, d2d3_pid);
		FillPairPid(event3, event4, config.track_.d4d3_window, d3d4_pid);
		FillSiliconPid(event4, event_s, d4s_pid);
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	d1d2_pid.Write();
	d2d3_pid.Write();
	d3d4_pid.Write();
	d4s_pid.Write();
	opf.Close();

	return 0;
}
