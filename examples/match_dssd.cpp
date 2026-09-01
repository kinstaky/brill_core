#include "include/t0/dssd.h"

#include <filesystem>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/dssd_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("match_dssd", "Match normalized DSSD events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("w,window", "Override matching tolerance.", cxxopts::value<double>(), "window")
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
		config.root.trigger = result["trigger"].as<std::string>();
	}
	const int run = result["run"].as<int>();
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}
	int normalize_file_run = 0;
	for (const auto &[start, use] : config.normalize.runs) {
		if (run >= start) normalize_file_run = use;
	}

	for (const std::string &detector_name : detectors) {
		const brill::SiliconDetectorConfig *detector =
			config.FindDetector(detector_name);
		if (!detector) {
			std::cerr << "Error: Detector " << detector_name << " is not found in config.\n";
			return 1;
		}
		double match_tolerance = detector->match_tolerance;
		if (result.count("window")) {
			match_tolerance = result["window"].as<double>();
		}
		brill::SiliconDetectorConfig working_detector = *detector;
		working_detector.match_tolerance = match_tolerance;

		brill::DssdNormalizeParameters parameters(detector->front_strips, detector->back_strips);
		std::string normalize_dir = brill::JoinPath(config.root.workspace, config.paths.normalize);
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
			brill::JoinPath(config.root.workspace, config.paths.ingot).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.root.trigger).c_str(),
			run
		);
		TFile ipf(input_path, "read");
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from " << input_path << " failed.\n";
			return 1;
		}
		brill::DssdEvent raw_event;
		brill::DssdEvent normalized_event;
		brill::SetupInput(ipt, raw_event);

		TString output_path = TString::Format(
			"%s/%s_%s%04d.root",
			brill::JoinPath(config.root.workspace, config.paths.match).c_str(),
			detector_name.c_str(),
			brill::TriggerInfix(config.root.trigger).c_str(),
			run
		);
		TFile opf(output_path, "recreate");
		TTree opt("tree", "matched dssd");
		brill::DssdMatchEvent match_event;
		brill::SetupOutput(&opt, match_event);

		long long total = ipt->GetEntries();
		long long last_percentage = 0;
		printf("Matching %s   0%%", detector_name.c_str());
		fflush(stdout);
		for (long long entry = 0; entry < total; ++entry) {
			if (entry * 100ll / total > last_percentage) {
				last_percentage = entry * 100ll / total;
				printf("\b\b\b\b%3lld%%", last_percentage);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			parameters.Apply(raw_event, normalized_event);
			brill::MatchDssdEvent(normalized_event, working_detector, match_event);
			opt.Fill();
		}
		printf("\b\b\b\b100%%\n");

		opf.cd();
		opt.Write();
		opf.Close();
		ipf.Close();
	}

	return 0;
}
