#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <string>

#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double AxisPitchX(const brill::SiliconDetectorConfig &detector) {
	if (detector.name == "t0d2") {
		return detector.size_x_mm / double(detector.back_strips);
	}
	return detector.size_x_mm / double(detector.front_strips);
}

double AxisPitchY(const brill::SiliconDetectorConfig &detector) {
	if (detector.name == "t0d2") {
		return detector.size_y_mm / double(detector.front_strips);
	}
	return detector.size_y_mm / double(detector.back_strips);
}

double RandomizeCoordinate(double value, double pitch, TRandom3 &generator) {
	return value + generator.Rndm() * pitch - 0.5 * pitch;
}

void FillOffsets(
	TH1F &hist_x,
	TH1F &hist_y,
	const brill::DssdMatchEvent &left,
	double left_pitch_x,
	double left_pitch_y,
	const brill::DssdMatchEvent &right,
	double right_pitch_x,
	double right_pitch_y,
	TRandom3 &generator
) {
	for (int i = 0; i < left.num; ++i) {
		double left_x = RandomizeCoordinate(left.x[i], left_pitch_x, generator);
		double left_y = RandomizeCoordinate(left.y[i], left_pitch_y, generator);
		for (int j = 0; j < right.num; ++j) {
			double right_x = RandomizeCoordinate(right.x[j], right_pitch_x, generator);
			double right_y = RandomizeCoordinate(right.y[j], right_pitch_y, generator);
			hist_x.Fill(right_x - left_x);
			hist_y.Fill(right_y - left_y);
		}
	}
}

void FillStripOffsets(
	TH1F &hist_x,
	TH1F &hist_y,
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right
) {
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			hist_x.Fill(right.x[i] - left.x[j]);
			hist_y.Fill(right.y[i] - left.y[j]);
		}
	}
}

void FitAndPrint(TH1F &histogram) {
	TF1 fit((std::string("fit_") + histogram.GetName()).c_str(), "gaus", -5.0, 5.0);
	fit.SetParameter(1, 0.0);
	fit.SetParameter(2, 1.0);
	if (histogram.GetEntries() > 0) {
		histogram.Fit(&fit, "QR+");
		std::cout
			<< histogram.GetName()
			<< ": mean = " << fit.GetParameter(1)
			<< ", sigma = " << fit.GetParameter(2)
			<< "\n";
	} else {
		std::cout
			<< histogram.GetName()
			<< ": empty histogram\n";
	}
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_center", "Estimate T0 detector center offsets.");
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
		config.root.trigger = result["trigger"].as<std::string>();
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const brill::SiliconDetectorConfig *detector1 = config.FindDetector("t0d1");
	const brill::SiliconDetectorConfig *detector2 = config.FindDetector("t0d2");
	const brill::SiliconDetectorConfig *detector3 = config.FindDetector("t0d3");
	const brill::SiliconDetectorConfig *detector4 = config.FindDetector("t0d4");
	if (!detector1 || !detector2 || !detector3 || !detector4) {
		std::cerr << "Error: Missing T0 detector config.\n";
		return 1;
	}

	const std::string match_dir = brill::JoinPath(
		config.root.workspace, config.paths.match
	);
	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);

	TChain chain1("tree");
	TChain chain2("tree");
	TChain chain3("tree");
	TChain chain4("tree");
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
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain3, "d3");
	chain1.AddFriend(&chain4, "d4");

	brill::DssdMatchEvent event1;
	brill::DssdMatchEvent event2;
	brill::DssdMatchEvent event3;
	brill::DssdMatchEvent event4;
	brill::SetupInput(&chain1, event1);
	brill::SetupInput(&chain1, event2, "d2.");
	brill::SetupInput(&chain1, event3, "d3.");
	brill::SetupInput(&chain1, event4, "d4.");

	TString output_path = TString::Format(
		"%s/t0_center_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TH1F d2d1_x_offset("d2d1x", "D2 x - D1 x", 1000, -10.0, 10.0);
	TH1F d2d1_y_offset("d2d1y", "D2 y - D1 y", 1000, -10.0, 10.0);
	TH1F d3d1_x_offset("d3d1x", "D3 x - D1 x", 1000, -10.0, 10.0);
	TH1F d3d1_y_offset("d3d1y", "D3 y - D1 y", 1000, -10.0, 10.0);
	TH1F d4d1_x_offset("d4d1x", "D4 x - D1 x", 1000, -10.0, 10.0);
	TH1F d4d1_y_offset("d4d1y", "D4 y - D1 y", 1000, -10.0, 10.0);
	TH1F d3d2_x_offset("d3d2x", "D3 x - D2 x", 1000, -10.0, 10.0);
	TH1F d3d2_y_offset("d3d2y", "D3 y - D2 y", 1000, -10.0, 10.0);
	TH1F d4d2_x_offset("d4d2x", "D4 x - D2 x", 1000, -10.0, 10.0);
	TH1F d4d2_y_offset("d4d2y", "D4 y - D2 y", 1000, -10.0, 10.0);
	TH1F d4d3_x_offset("d4d3x", "D4 x - D3 x", 1000, -10.0, 10.0);
	TH1F d4d3_y_offset("d4d3y", "D4 y - D3 y", 1000, -10.0, 10.0);
	TH1F d2d1xs("d2d1xs", "D2 x - D1 x", 40, -10.0, 10.0);
	TH1F d2d1ys("d2d1ys", "D2 y - D1 y", 40, -10.0, 10.0);
	TH1F d3d1xs("d3d1xs", "D3 x - D1 x", 40, -10.0, 10.0);
	TH1F d3d1ys("d3d1ys", "D3 y - D1 y", 40, -10.0, 10.0);
	TH1F d4d1xs("d4d1xs", "D4 x - D1 x", 40, -10.0, 10.0);
	TH1F d4d1ys("d4d1ys", "D4 y - D1 y", 40, -10.0, 10.0);
	TH1F d3d2xs("d3d2xs", "D3 x - D2 x", 40, -10.0, 10.0);
	TH1F d3d2ys("d3d2ys", "D3 y - D2 y", 40, -10.0, 10.0);
	TH1F d4d2xs("d4d2xs", "D4 x - D2 x", 40, -10.0, 10.0);
	TH1F d4d2ys("d4d2ys", "D4 y - D2 y", 40, -10.0, 10.0);
	TH1F d4d3xs("d4d3xs", "D4 x - D3 x", 40, -10.0, 10.0);
	TH1F d4d3ys("d4d3ys", "D4 y - D3 y", 40, -10.0, 10.0);

	const double d1_pitch_x = AxisPitchX(*detector1);
	const double d1_pitch_y = AxisPitchY(*detector1);
	const double d2_pitch_x = AxisPitchX(*detector2);
	const double d2_pitch_y = AxisPitchY(*detector2);
	const double d3_pitch_x = AxisPitchX(*detector3);
	const double d3_pitch_y = AxisPitchY(*detector3);
	const double d4_pitch_x = AxisPitchX(*detector4);
	const double d4_pitch_y = AxisPitchY(*detector4);

	const long long total = chain1.GetEntries();
	TRandom3 generator(total);
	long long last_percentage = -1;
	std::printf("Estimating T0 center   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);

		FillOffsets(
			d2d1_x_offset,
			d2d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			generator
		);
		FillOffsets(
			d3d1_x_offset,
			d3d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			generator
		);
		FillOffsets(
			d4d1_x_offset,
			d4d1_y_offset,
			event1,
			d1_pitch_x,
			d1_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
		FillOffsets(
			d3d2_x_offset,
			d3d2_y_offset,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			generator
		);
		FillOffsets(
			d4d2_x_offset,
			d4d2_y_offset,
			event2,
			d2_pitch_x,
			d2_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
		FillOffsets(
			d4d3_x_offset,
			d4d3_y_offset,
			event3,
			d3_pitch_x,
			d3_pitch_y,
			event4,
			d4_pitch_x,
			d4_pitch_y,
			generator
		);
		FillStripOffsets(d2d1xs, d2d1ys, event1, event2);
		FillStripOffsets(d3d1xs, d3d1ys, event1, event3);
		FillStripOffsets(d4d1xs, d4d1ys, event1, event4);
		FillStripOffsets(d3d2xs, d3d2ys, event2, event3);
		FillStripOffsets(d4d2xs, d4d2ys, event2, event4);
		FillStripOffsets(d4d3xs, d4d3ys, event3, event4);
	}
	std::printf("\b\b\b\b100%%\n");

	FitAndPrint(d2d1_x_offset);
	FitAndPrint(d2d1_y_offset);
	FitAndPrint(d3d1_x_offset);
	FitAndPrint(d3d1_y_offset);
	FitAndPrint(d4d1_x_offset);
	FitAndPrint(d4d1_y_offset);
	FitAndPrint(d3d2_x_offset);
	FitAndPrint(d3d2_y_offset);
	FitAndPrint(d4d2_x_offset);
	FitAndPrint(d4d2_y_offset);
	FitAndPrint(d4d3_x_offset);
	FitAndPrint(d4d3_y_offset);

	opf.cd();
	d2d1_x_offset.Write();
	d2d1_y_offset.Write();
	d3d1_x_offset.Write();
	d3d1_y_offset.Write();
	d4d1_x_offset.Write();
	d4d1_y_offset.Write();
	d3d2_x_offset.Write();
	d3d2_y_offset.Write();
	d4d2_x_offset.Write();
	d4d2_y_offset.Write();
	d4d3_x_offset.Write();
	d4d3_y_offset.Write();
	d2d1xs.Write();
	d2d1ys.Write();
	d3d1xs.Write();
	d3d1ys.Write();
	d4d1xs.Write();
	d4d1ys.Write();
	d3d2xs.Write();
	d3d2ys.Write();
	d4d2xs.Write();
	d4d2ys.Write();
	d4d3xs.Write();
	d4d3ys.Write();
	opf.Close();

	return 0;
}
