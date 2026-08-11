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
#include "include/event/ingot/csi_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

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

	const std::string ingot_dir = brill::JoinPath(config.root_.workspace, config.paths_.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.root_.trigger);

	TChain chain_s("tree");
	TChain chain_csi("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
		++added_runs;
		chain_s.Add(TString::Format(
			"%s/t0s_%s%04d.root",
			ingot_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain_csi.Add(TString::Format(
			"%s/t0csi_%s%04d.root",
			ingot_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}
	chain_s.AddFriend(&chain_csi, "csi");

	brill::SiliconEvent ssd;
	brill::CsiEvent<36> csi;
	brill::SetupInput(&chain_s, ssd);
	brill::SetupInput(&chain_s, csi, "csi.");

	TString output_path = TString::Format(
		"%s/t0csi_pid_%s%04d_%04d.root",
		brill::JoinPath(config.root_.workspace, config.paths_.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");

	TH2F* pids[36];
	double csi_range[36] = {
		5000.0, 5000.0, 5000.0, 5000.0,
		5000.0, 5000.0, 5000.0, 5000.0,
		5000.0, 5000.0, 5000.0, 5000.0,
		5000.0, 5000.0, 5000.0, 12000.0,
		12000.0, 12000.0, 12000.0, 12000.0,
		32000.0, 32000.0, 32000.0, 64000.0,
		32000.0, 32000.0, 32000.0, 32000.0,
		32000.0, 12000.0, 32000.0, 32000.0,
		32000.0, 32000.0, 32000.0, 32000.0
	};
	for (int i = 0; i < 36; ++i) {
		pids[i] = new TH2F(
			TString::Format("h%d", i),
			TString::Format("CsI %d pid", i),
			500.0, 0.0, csi_range[i], 500.0, 0.0, 30000.0
		);
	}

	const long long total = chain_s.GetEntries();
	long long last_percentage = 0;
	std::printf("Filling T0 CsI PID   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain_s.GetEntry(entry);
		if (!ssd.valid) continue;
		for (int i = 0; i < 36; ++i) {
			if (!csi.valid[i]) continue;
			pids[i]->Fill(csi.energy[i], ssd.energy);
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	for (int i = 0; i < 36; ++i) pids[i]->Write();
	opf.Close();

	return 0;
}
