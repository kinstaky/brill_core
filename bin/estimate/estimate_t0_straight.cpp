#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"

namespace {

struct ParticleStraightSpec {
	std::string particle;
	bool tail = false;
	std::unique_ptr<TCutG> cut;
	TGraph graph;
	std::unique_ptr<TF1> fit;
	double a = 0.5;
	double b = -0.04;
	double c = 0.0;
};

struct SliceStraightSpec {
	std::string key;
	const brill::TrackWindowConfig *window = nullptr;
	bool use_t0s = false;
	std::vector<ParticleStraightSpec> particles;
	std::unique_ptr<TH2F> pid;
	std::unique_ptr<TH2F> straight_pid;
	std::unique_ptr<TH1F> ef_hist;
	double default_a = 0.5;
	double default_b = -0.04;
};

ParticleStraightSpec MakeParticle(const std::string &particle, bool tail = false) {
	ParticleStraightSpec spec;
	spec.particle = particle;
	spec.tail = tail;
	return spec;
}

SliceStraightSpec MakeSlice(
	const std::string &key,
	const brill::TrackWindowConfig *window,
	bool use_t0s,
	std::unique_ptr<TH2F> pid,
	std::unique_ptr<TH2F> straight_pid,
	std::unique_ptr<TH1F> ef_hist
) {
	SliceStraightSpec spec;
	spec.key = key;
	spec.window = window;
	spec.use_t0s = use_t0s;
	spec.pid = std::move(pid);
	spec.straight_pid = std::move(straight_pid);
	spec.ef_hist = std::move(ef_hist);
	return spec;
}


void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double PidFit(double *x, double *par) {
	return (0.5 / par[0]) * (
		std::sqrt(std::pow(x[0], 2.0) + 4.0 * par[0] * std::pow(par[1] * x[0] - par[2], 2.0))
		- x[0]
	);
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

void FillGraphsForPoint(SliceStraightSpec &slice, double e, double de) {
	for (auto &particle : slice.particles) {
		if (particle.cut && particle.cut->IsInside(e, de)) {
			particle.graph.AddPoint(e, de);
		}
	}
}

void FillStraightPoint(SliceStraightSpec &slice, double e, double de) {
	double ef = std::sqrt(de * e + slice.default_a * de * de) + slice.default_b * e;
	slice.pid->Fill(e, de);
	slice.straight_pid->Fill(e, ef);
	slice.ef_hist->Fill(ef);
}

template <typename Handler>
void HandleDssdPairs(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	const brill::TrackWindowConfig &window,
	Handler handler
) {
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			if (!InTrackWindow(left, i, right, j, window)) continue;
			handler(right.energy[j], left.energy[i]);
		}
	}
}

template <typename Handler>
void HandleSiliconPairs(
	const brill::DssdMatchEvent &dssd,
	const brill::SiliconEvent &silicon,
	Handler handler
) {
	if (!silicon.valid) return;
	for (int i = 0; i < dssd.num; ++i) {
		handler(double(silicon.energy), dssd.energy[i]);
	}
}

void FitSlice(SliceStraightSpec &slice) {
	std::cout << slice.key << " straight PID parameters.\n";
	for (auto &particle : slice.particles) {
		std::string graph_name = "g" + slice.key + "_" + particle.particle + (particle.tail ? "_tail" : "");
		particle.graph.SetName(graph_name.c_str());
		particle.graph.SetTitle(graph_name.c_str());
		if (particle.graph.GetN() == 0 || !particle.cut) {
			std::cout << "  " << particle.particle << ": empty\n";
			continue;
		}

		double x_min = particle.cut->GetX()[0];
		double x_max = particle.cut->GetX()[0];
		double y_max = particle.cut->GetY()[0];
		for (int i = 1; i < particle.cut->GetN(); ++i) {
			if (particle.cut->GetX()[i] < x_min) x_min = particle.cut->GetX()[i];
			if (particle.cut->GetX()[i] > x_max) x_max = particle.cut->GetX()[i];
			if (particle.cut->GetY()[i] > y_max) y_max = particle.cut->GetY()[i];
		}

		std::string fit_name = "f" + slice.key + "_" + particle.particle + (particle.tail ? "_tail" : "");
		particle.fit = std::make_unique<TF1>(fit_name.c_str(), PidFit, x_min, x_max, 3);
		particle.fit->SetParameter(0, particle.a);
		particle.fit->SetParameter(1, particle.b);
		particle.fit->SetParameter(2, y_max);
		particle.fit->SetParLimits(2, 0.0, 1e10);
		particle.graph.Fit(particle.fit.get(), "RQ+");
		particle.a = particle.fit->GetParameter(0);
		particle.b = particle.fit->GetParameter(1);
		particle.c = particle.fit->GetParameter(2);
		std::cout
			<< "  " << particle.particle
			<< ": A " << particle.a
			<< ", B " << particle.b
			<< ", C " << particle.c
			<< "\n";

		if (particle.particle == "4He" && !particle.tail) {
			slice.default_a = particle.a;
			slice.default_b = particle.b;
		}
	}
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_straight", "Estimate T0 straight PID.");
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

	std::vector<SliceStraightSpec> slices;
	slices.push_back(MakeSlice(
		"t0d1d2",
		&config.track_.d2d1_window,
		false,
		std::make_unique<TH2F>("d1d2", "D1-D2 PID", 5000, 0.0, 60000.0, 5000, 0.0, 30000.0),
		std::make_unique<TH2F>("sd1d2", "D1-D2 straight PID", 5000, 0.0, 60000.0, 5000, 0.0, 60000.0),
		std::make_unique<TH1F>("ed1d2", "D1-D2 straight energy", 5000, 0.0, 60000.0)
	));
	slices.back().particles.push_back(MakeParticle("4He"));

	slices.push_back(MakeSlice(
		"t0d2d3",
		&config.track_.d3d2_window,
		false,
		std::make_unique<TH2F>("d2d3", "D2-D3 PID", 5000, 0.0, 50000.0, 5000, 0.0, 80000.0),
		std::make_unique<TH2F>("sd2d3", "D2-D3 straight PID", 5000, 0.0, 50000.0, 5000, 0.0, 50000.0),
		std::make_unique<TH1F>("ed2d3", "D2-D3 straight energy", 5000, 0.0, 50000.0)
	));
	for (const char *particle : {"4He", "6Li", "7Be", "10B", "10C"}) {
		slices.back().particles.push_back(MakeParticle(particle));
	}

	slices.push_back(MakeSlice(
		"t0d3d4",
		&config.track_.d4d3_window,
		false,
		std::make_unique<TH2F>("d3d4", "D3-D4 PID", 5000, 0.0, 45000.0, 5000, 0.0, 45000.0),
		std::make_unique<TH2F>("sd3d4", "D3-D4 straight PID", 5000, 0.0, 45000.0, 5000, 0.0, 45000.0),
		std::make_unique<TH1F>("ed3d4", "D3-D4 straight energy", 5000, 0.0, 45000.0)
	));
	for (const char *particle : {"4He", "6Li", "7Be", "10B", "10C"}) {
		slices.back().particles.push_back(MakeParticle(particle));
	}

	slices.push_back(MakeSlice(
		"t0d4s",
		nullptr,
		true,
		std::make_unique<TH2F>("d4s", "D4-S PID", 5000, 0.0, 80000.0, 5000, 0.0, 45000.0),
		std::make_unique<TH2F>("sd4s", "D4-S straight PID", 5000, 0.0, 80000.0, 5000, 0.0, 80000.0),
		std::make_unique<TH1F>("ed4s", "D4-S straight energy", 5000, 0.0, 80000.0)
	));
	for (const char *particle : {"4He", "6Li"}) {
		slices.back().particles.push_back(MakeParticle(particle));
	}

	for (auto &slice : slices) {
		for (auto &particle : slice.particles) {
			if (brill::ParseCutFile(
				config.root_.workspace,
				slice.key,
				particle.particle,
				particle.tail,
				particle.cut
			)) {
				return 1;
			}
		}
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
		"%s/t0_straight_%s%04d_%04d.root",
		brill::JoinPath(config.root_.workspace, config.paths_.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Filling T0 straight fit points   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		HandleDssdPairs(event1, event2, config.track_.d2d1_window, [&](double e, double de) {
			FillGraphsForPoint(slices[0], e, de);
		});
		HandleDssdPairs(event2, event3, config.track_.d3d2_window, [&](double e, double de) {
			FillGraphsForPoint(slices[1], e, de);
		});
		HandleDssdPairs(event3, event4, config.track_.d4d3_window, [&](double e, double de) {
			FillGraphsForPoint(slices[2], e, de);
		});
		HandleSiliconPairs(event4, event_s, [&](double e, double de) {
			FillGraphsForPoint(slices[3], e, de);
		});
	}
	std::printf("\b\b\b\b100%%\n");

	for (auto &slice : slices) {
		FitSlice(slice);
	}

	slices[0].default_a = 0.1128;
	slices[0].default_b = -0.0738;
	slices[1].default_a = 0.4656;
	slices[1].default_b = -0.0513;
	slices[2].default_a = 0.7022;
	slices[2].default_b = -0.0491;
	slices[3].default_a = 0.9458;
	slices[3].default_b = -0.0429;

	last_percentage = -1;
	std::printf("Filling T0 straight PID   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		HandleDssdPairs(event1, event2, config.track_.d2d1_window, [&](double e, double de) {
			FillStraightPoint(slices[0], e, de);
		});
		HandleDssdPairs(event2, event3, config.track_.d3d2_window, [&](double e, double de) {
			FillStraightPoint(slices[1], e, de);
		});
		HandleDssdPairs(event3, event4, config.track_.d4d3_window, [&](double e, double de) {
			FillStraightPoint(slices[2], e, de);
		});
		HandleSiliconPairs(event4, event_s, [&](double e, double de) {
			FillStraightPoint(slices[3], e, de);
		});
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	for (auto &slice : slices) {
		for (auto &particle : slice.particles) {
			if (particle.cut) particle.cut->Write();
		}
	}
	for (auto &slice : slices) {
		for (auto &particle : slice.particles) {
			particle.graph.Write();
		}
	}
	for (auto &slice : slices) slice.pid->Write();
	for (auto &slice : slices) slice.straight_pid->Write();
	for (auto &slice : slices) slice.ef_hist->Write();
	opf.Close();

	return 0;
}
