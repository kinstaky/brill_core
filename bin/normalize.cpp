#include "include/t0/dssd.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TGraph.h>
#include <TString.h>
#include <TF1.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/t0/dssd.h"
#include "include/utils.h"

int NormalizeStrips(
	const brill::NromalizeStripsConfig &config,
	const bool use_integral,
	TChain &chain,
	const brill::DssdEvent &event,
	brill::DssdNormalizeParameters &parameters,
	std::vector<bool> &has_normalized
) {
	const int &side = config.norm_side;
	const int offset = side * parameters.front_strips_;
	// energy graph fe:be or be:fe
	TGraph ge[128];

	// total number of entries
	long long total = chain.GetEntries();
	// 1/100 of total entries
	long long last_percentage = 0;
	// show start
	printf(
		"Filling for side %d [%d, %d] reference strips [%d, %d]   0%%",
		config.norm_side, config.norm[0], config.norm[1], config.ref[0], config.ref[1]
	);
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		// show process
		if (entry * 100ll / total > last_percentage) {
			last_percentage = entry * 100ll / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		chain.GetEntry(entry);
		// ignore multiple hit events
		if (event.front_num != 1 || event.back_num != 1) continue;

		const int &fs = event.front_strip[0];
		const int &bs = event.back_strip[0];
		const double &fe = use_integral ? event.front_integral[0] : event.front_energy[0];
		const double &be = use_integral ? event.back_integral[0] : event.back_energy[0];

		if (config.norm_side == 0) {
			// jump if not reference strips
			if (bs < config.ref[0] || bs > config.ref[1]) continue;
			// jump if not normalize strips
			if (fs < config.norm[0] || fs > config.norm[1]) continue;
			// jump if has normalized
			if (has_normalized[fs]) continue;
			// jump if energy out of range
			if (be < config.ref_energy[0] || be > config.ref_energy[1]) continue;
			if (fe < config.norm_energy[0] || fe > config.norm_energy[1]) continue;
			// fill to graph
			ge[fs].AddPoint(fe, parameters.NormEnergy(1, bs, be));
		} else {
			// jump if not reference strips
			if (fs < config.ref[0]|| fs > config.ref[1]) continue;
			// jump if not normalize strips
			if (bs < config.norm[0] || bs > config.norm[1]) continue;
			// jump if has normalized
			if (has_normalized[offset+bs]) continue;
			// jump if energy out of range
			if (fe < config.ref_energy[0] || fe > config.ref_energy[1]) continue;
			if (be < config.norm_energy[0] || be > config.norm_energy[1]) continue;
	/*std::cout << config.ref[0] << ", " << config.ref[1] << ", "
		<< config.norm[0] << ", " << config.norm[1] << ", "
		<< config.ref_energy[0] << ", " << config.ref_energy[1] << ", "
		<< config.norm_energy[0] << ", " << config.norm_energy[1] << ", "
		<< fe << ", " << be << ", " << fs << ", " << bs << "\n";*/
			// fill to graph
			ge[bs].AddPoint(be, parameters.NormEnergy(0, fs, fe));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fitting
	std::cout << "side " << side << " normalize parameters.\n";
	for (int i = config.norm[0]; i <= config.norm[1]; ++i) {
		if (has_normalized[offset+i]) continue;
		// only fits when over 10 points
		if (ge[i].GetN() > 5) {
			// fitting function
			TF1 energy_fit("efit", "pol1", 0, 60000);
			// set initial value
			energy_fit.SetParameter(0, 0.0);
			energy_fit.SetParameter(1, 1.0);
			//energy_fit.SetParLimits(2, -1e-6, 1e-6);
			//energy_fit.SetParameter(2, 0.0);
			// fit
			ge[i].Fit(&energy_fit, "QR+ ROB=0.8");
			// store the normalized parameters
			if (config.norm_side == 0) {
				parameters.front_p0_[i] = energy_fit.GetParameter(0);
				parameters.front_p1_[i] = energy_fit.GetParameter(1);
				//parameters.front_p2[i] = energy_fit.GetParameter(2);
			} else {
				parameters.back_p0_[i] = energy_fit.GetParameter(0);
				parameters.back_p1_[i] = energy_fit.GetParameter(1);
				//parameters.back_p2[i] = energy_fit.GetParameter(2);
			}
		}
		// store the graph
		ge[i].Write(TString::Format("g%c%d", "fb"[side], i));
		// set as normalized
		has_normalized[offset+i] = true;
		// print normalized paramters on screen
		if (side == 0) {
			std::cout << i
				<< " " << parameters.front_p0_[i]
				<< ", " << parameters.front_p1_[i]
				<< ", " << parameters.front_p2_[i]
				<< "\n";
		} else {
			std::cout << i
				<< " " << parameters.back_p0_[i]
				<< ", " << parameters.back_p1_[i]
				<< ", " << parameters.back_p2_[i]
				<< "\n";
		}
	}

	// residual
	TGraph res[128];
	for (int i = config.norm[0]; i < config.norm[1]; ++i) {
		if (has_normalized[offset+i]) continue;
		int point = ge[i].GetN();
		double *gex = ge[i].GetX();
		double *gey = ge[i].GetY();
		for (int j = 0; j < point; ++j) {
			res[i].AddPoint(gex[j], parameters.NormEnergy(side, i, gex[j])-gey[j]);
		}
		res[i].Write(TString::Format("res%c%d", "fb"[side], i));
	}
	return 0;

}
void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("match_dssd", "Match normalized DSSD events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"detector",
			"Detectors to match.",
			cxxopts::value<std::vector<std::string>>(),
			"detector"
		);
	options.parse_positional({"detector"});

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
	if (!result.count("detector")) {
		std::cerr << "Error: Missing detector positional arguments.\n";
		PrintUsage(options);
		return 1;
	}

	const std::set<std::string> allowed_detectors = {"t0d1", "t0d2", "t0d3", "t0d4"};
	std::vector<std::string> detectors = result["detector"].as<std::vector<std::string>>();
	for (const std::string &detector : detectors) {
		if (allowed_detectors.count(detector) == 0) {
			std::cerr << "Error: Unsupported detector " << detector << ".\n";
			return 1;
		}
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
		return -1;
	}

	for (const std::string &detector_name : detectors) {
		const brill::SiliconDetectorConfig *detector =
			config.FindDetector(detector_name);
		if (!detector) {
			std::cerr << "Error: Detector " << detector_name << " is not found in config.\n";
			return 1;
		}

		TChain chain("tree");
		int added_runs = 0;
		for (int current_run = run; current_run <= end_run; ++current_run) {
			if (config.IsJumpRun(current_run)) continue;
			++added_runs;
			chain.Add(TString::Format(
				"%s/%s_%s%04d.root",
				brill::JoinPath(config.root_.workspace, config.paths_.ingot).c_str(),
				detector_name.c_str(),
				brill::TriggerInfix(config.root_.trigger).c_str(),
				current_run
			));
		}
		if (added_runs == 0) {
			std::cout << "No runs to process after jumping runs.\n";
			return 0;
		}

		brill::DssdNormalizeParameters parameters(detector->front_strips, detector->back_strips);

		brill::DssdEvent raw_event;
		brill::SetupInput(&chain, raw_event);

		std::string normalize_dir = brill::JoinPath(config.root_.workspace, config.paths_.normalize);
		TString output_path = TString::Format(
			"%s/%s_%s%04d_%04d.root",
			normalize_dir.c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.root_.trigger).c_str(),
			run,
			end_run
		);
		TFile opf(output_path, "recreate");

		const auto &strips_config = config.normalize_.detectors[detector_name].strips;

		std::vector<bool> has_normalized;
		for (int i = 0; i < parameters.front_strips_+parameters.back_strips_; ++i) {
			has_normalized.push_back(false);
		}

		opf.cd();
		for (const auto &strips : strips_config) {
			NormalizeStrips(
				strips, detector->use_integral, chain, raw_event,
				parameters, has_normalized
			);
		}

		TString parameter_path = TString::Format(
			"%s/%s_%04d.txt",
			normalize_dir.c_str(),
			detector_name.c_str(),
			run
		);
		if (parameters.Write(parameter_path.Data())) {
			std::cerr << "Error: Write normalize parameters failed.\n";
			return 1;
		}

		opf.cd();
		opf.Close();
	}

	return 0;
}
