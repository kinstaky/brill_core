#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/lost_energy_calculator.h"
#include "include/event/particle_event.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"

namespace {

constexpr int kCalibrationLayers = 5;
constexpr int kFirstStopLayer = 2;
constexpr int kLastStopLayer = 5;
constexpr int kD2Bit = 0x2;
constexpr int kD3Bit = 0x4;
constexpr int kD4Bit = 0x8;
constexpr int kSBit = 0x10;

struct CalibrationParameters {
	double p0[kCalibrationLayers] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double p1[kCalibrationLayers] = {1.0, 1.0, 1.0, 1.0, 1.0};
};

struct ParticleEntry {
	int charge = 0;
	int mass = 0;
	double energy = 0.0;
	double time = 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	bool stop = false;
	int last = -1;
};

using StopEnergyCache = std::map<std::tuple<int, int, int>, double>;

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int ReadCalibrationParameters(
	const std::string &path,
	CalibrationParameters &parameters
) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open calibration parameter file " << path << " failed.\n";
		return -1;
	}

	std::string header;
	std::getline(fin, header);
	int index = -1;
	double p0 = 0.0;
	double p1 = 1.0;
	while (fin >> index >> p0 >> p1) {
		if (index < 0 || index >= kCalibrationLayers) continue;
		parameters.p0[index] = p0;
		parameters.p1[index] = p1;
	}
	return 0;
}

double CalibrateEnergy(
	const CalibrationParameters &parameters,
	int layer,
	double raw_energy
) {
	return parameters.p0[layer] + parameters.p1[layer] * raw_energy;
}

double RebuildEnergy(
	const brill::T0Event &input,
	int particle_index,
	const CalibrationParameters &parameters,
	const double max_stop_energy
) {
	double energy = 0.0;
	int actual_layer = input.layer[particle_index] == 6 ? 5 : input.layer[particle_index];
	for (int i = 1; i < actual_layer; ++i) {
		double layer_energy = CalibrateEnergy(parameters, i, input.energy[particle_index][i]);
		if (i == input.layer[particle_index] - 1 && layer_energy > max_stop_energy) {
			return 0.0;
		}
		energy += layer_energy;
	}
	return energy;
}

std::string RangeCachePath(
	const brill::AppConfig &config,
	int charge,
	int mass
) {
	return TString::Format(
		"%s/si_z%d_a%d.root",
		brill::JoinPath(config.root.workspace, config.paths.energy_calculator).c_str(),
		charge,
		mass
	).Data();
}


double MaximumStopEnergy(
	const brill::AppConfig &config,
	int charge,
	int mass,
	int stop_layer,
	StopEnergyCache &cache
) {
	if (stop_layer > 5) return 0.0;
	const auto key = std::make_tuple(charge, mass, stop_layer);
	auto iterator = cache.find(key);
	if (iterator != cache.end()) return iterator->second;

	if (
		stop_layer < kFirstStopLayer
		|| stop_layer > kLastStopLayer
		|| size_t(stop_layer - 1) >= config.t0.silicon.size()
	) {
		return 0.0;
	}

	const std::string &detector_name = config.t0.silicon[stop_layer - 1];
	const auto *detector = config.FindDetector(detector_name);
	if (!detector) {
		std::cerr
			<< "Error: Detector config for "
			<< detector_name
			<< " not found when rebuilding T0.\n";
		return 0.0;
	}

	brill::LostEnergyCalculator calculator(
		charge,
		mass,
		brill::SiliconMaterial(),
		detector->thickness_um,
		RangeCachePath(config, charge, mass)
	);
	double maximum_stop_energy = calculator.IncidentEnergy(0.0);
	cache[key] = maximum_stop_energy;
	return maximum_stop_energy;
}


bool CompareParticleEntry(const ParticleEntry &left, const ParticleEntry &right) {
	if (left.charge != right.charge) return left.charge < right.charge;
	if (left.mass != right.mass) return left.mass < right.mass;
	return left.last < right.last;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("rebuild_t0", "Rebuild T0 particles from tracked events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
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
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);
	const std::string tele_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.telescope).c_str(),
		trigger_infix.c_str(),
		run
	).Data();
	const std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.root.workspace, config.paths.calibration).c_str()
	).Data();
	const std::string output_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.particle).c_str(),
		trigger_infix.c_str(),
		run
	).Data();

	CalibrationParameters calibration;
	if (ReadCalibrationParameters(calibration_path, calibration)) {
		return 1;
	}

	TFile ipf(tele_path.c_str(), "read");
	if (ipf.IsZombie()) {
		std::cerr << "Error: Open tracked T0 file " << tele_path << " failed.\n";
		return 1;
	}
	TTree *ipt = static_cast<TTree*>(ipf.Get("tree"));
	if (!ipt) {
		std::cerr << "Error: Get tree from " << tele_path << " failed.\n";
		return 1;
	}

	std::filesystem::path output_file_path(output_path);
	if (!output_file_path.parent_path().empty()) {
		std::filesystem::create_directories(output_file_path.parent_path());
	}
	TFile opf(output_path.c_str(), "recreate");
	TTree opt("tree", "rebuild t0 particles");

	brill::T0Event input_event;
	brill::ParticleEvent output_event;
	brill::SetupInput(ipt, input_event);
	brill::Reset(output_event);
	brill::SetupOutput(&opt, output_event);
	StopEnergyCache stop_energy_cache;

	const long long total = ipt->GetEntries();
	long long last_percentage = -1;
	std::printf("Rebuilding T0   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}

		ipt->GetEntry(entry);
		brill::Reset(output_event);
		std::vector<ParticleEntry> particles;
		for (int i = 0; i < input_event.num; ++i) {
			if (input_event.charge[i] <= 0 || input_event.mass[i] <= 0) continue;
			double max_stop_energy = MaximumStopEnergy(
				config,
				input_event.charge[i],
				input_event.mass[i],
				input_event.layer[i],
				stop_energy_cache
			);
			double particle_energy = RebuildEnergy(input_event, i, calibration, max_stop_energy);
			if (particle_energy == 0.0) continue;
			ParticleEntry particle;
			particle.charge = input_event.charge[i];
			particle.mass = input_event.mass[i];
			particle.energy = particle_energy;
			particle.stop = input_event.layer[i] != 6;
			particle.last = i;
			particle.x = input_event.x[i][1];
			particle.y = input_event.y[i][1];
			particle.z = input_event.z[i][1];
			particle.time = input_event.time[i][1];
			particles.push_back(particle);
		}

		std::sort(particles.begin(), particles.end(), CompareParticleEntry);
		for (const auto &particle : particles) {
			if (output_event.num >= 8) break;
			int index = output_event.num++;
			output_event.charge[index] = particle.charge;
			output_event.mass[index] = particle.mass;
			output_event.energy[index] = particle.energy;
			output_event.time[index] = particle.time;
			output_event.x[index] = particle.x;
			output_event.y[index] = particle.y;
			output_event.z[index] = particle.z;
			output_event.stop[index] = particle.stop;
			output_event.last[index] = particle.last;
		}

		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	ipf.Close();
	return 0;
}
