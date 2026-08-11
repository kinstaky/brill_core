#include <filesystem>
#include <cstdio>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/t0/dssd.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

void RemoveAdjacentStrips(brill::DssdEvent &event) {
	std::vector<int> front_index, back_index;
	for (int i = 0; i < event.front_num; ++i) {
		bool found = false;
		for (int j = i+1; j < event.front_num; ++j) {
			if (abs(event.front_strip[i]-event.front_strip[j])==1) {
				found = true;
				break;
			}
		}
		if (!found) front_index.push_back(i);
	}
	for (int i = 0; i < event.back_num; ++i) {
		bool found = false;
		for (int j = i+1; j < event.back_num; ++j) {
			if (abs(event.back_strip[i]-event.back_strip[j])==1) {
				found = true;
				break;
			}
		}
		if (!found) back_index.push_back(i);
	}
	event.front_num = int(front_index.size());
	event.back_num = int(back_index.size());
	for (size_t i = 0; i < front_index.size(); ++i) {
		event.front_strip[i] = event.front_strip[front_index[i]];
		event.front_energy[i] = event.front_energy[front_index[i]];
		event.front_time[i] = event.front_time[front_index[i]];
	}
	for (size_t i = 0; i < back_index.size(); ++i) {
		event.back_strip[i] = event.back_strip[back_index[i]];
		event.back_energy[i] = event.back_energy[back_index[i]];
		event.back_time[i] = event.back_time[back_index[i]];
	}
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_normalize", "Estimate normalized DSSD front-back result.");
	options.add_options()
		(
			"h,help",
			"Print help information."
		)
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"r,run",
			"Run number.",
			cxxopts::value<int>(),
			"run"
		)
		(
			"t,trigger",
			"Trigger type.",
			cxxopts::value<std::string>(),
			"trigger"
		)
		(
			"detector",
			"Detectors to estimate.",
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
	const int run = result["run"].as<int>();
	// const int end_run = result["end-run"].as<int>();
	// if (end_run < run) {
	// 	std::cout << "Error: end_run " << end_run << " should larger than run " << run << ".\n";
	// 	PrintUsage(options);
	// 	return 0;
	// }
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}
	if (result.count("trigger")) {
		config.root_.trigger = result["trigger"].as<std::string>();
	}
	int normalize_file_run = 0;
	for (const auto &[start, use] : config.normalize_.runs) {
		if (run >= start) normalize_file_run = use;
	}

	for (const std::string &detector_name : detectors) {
		const brill::SiliconDetectorConfig *detector =
			config.FindDetector(detector_name);
		if (!detector) {
			std::cerr << "Error: Detector " << detector_name << " is not found in config.\n";
			return 1;
		}

		brill::DssdNormalizeParameters parameters(detector->front_strips, detector->back_strips);

		std::string normalize_dir = brill::JoinPath(config.root_.workspace, config.paths_.normalize);
		TString parameter_path = TString::Format(
			"%s/%s_%04d.txt",
			normalize_dir.c_str(),
			detector_name.c_str(),
			normalize_file_run
		);
		if (parameters.Read(parameter_path.Data())) {
			return 1;
		}

		TString input_path = TString::Format(
			"%s/%s_%s%04d.root",
			brill::JoinPath(config.root_.workspace, config.paths_.ingot).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.root_.trigger).c_str(),
			run
		);
		TFile input_file(input_path, "read");
		TTree *input_tree = (TTree*)input_file.Get("tree");
		if (!input_tree) {
			std::cerr << "Error: Get tree from " << input_path << " failed.\n";
			return 1;
		}
		brill::DssdEvent raw_event;
		brill::SetupInput(input_tree, raw_event);

		TString output_path = TString::Format(
			"%s/%s_normalize_%s%04d.root",
			brill::JoinPath(config.root_.workspace, config.paths_.estimate).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.root_.trigger).c_str(),
			run
		);
		TFile output_file(output_path, "recreate");
		TTree opt("tree", "normalized dssd");
		brill::DssdEvent normalized_event;
		SetupOutput(&opt, normalized_event);
		TH1F hist("hde", "normalized front-back energy difference", 1000, -5000, 5000);

		printf("Estimating %s   0%%", detector_name.c_str());
		fflush(stdout);
		long long total = input_tree->GetEntries();
		long long last_percentage = 0;
		for (long long entry = 0; entry < total; ++entry) {
			if (entry * 100ll / total > last_percentage) {
				last_percentage = entry * 100ll / total;
				printf("\b\b\b\b%3lld%%", last_percentage);
				fflush(stdout);
			}
			input_tree->GetEntry(entry);
			if (detector->use_integral) {
				for (int i = 0; i < raw_event.front_num; ++i) {
					raw_event.front_energy[i] = raw_event.front_integral[i];
				}
				for (int i = 0; i < raw_event.back_num; ++i) {
					raw_event.back_energy[i] = raw_event.back_integral[i];
				}
			}
			parameters.Apply(raw_event, normalized_event);
			opt.Fill();
			RemoveAdjacentStrips(normalized_event);
			for (
				int i = 0;
				i < normalized_event.front_num && i < normalized_event.back_num;
				++i
			) {
				if (normalized_event.front_energy[i] == 0 || normalized_event.back_energy[i] == 0) continue;
				hist.Fill(normalized_event.front_energy[i] - normalized_event.back_energy[i]);
			}
		}
		printf("\b\b\b\b100%%\n");

		output_file.cd();
		hist.Write();
		opt.Write();
		output_file.Close();
		input_file.Close();
	}

	return 0;
}
