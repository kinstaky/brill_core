#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <Math/Vector3D.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/lost_energy_calculator.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/rebuild/nuclear_data.h"
#include "include/rebuild/particle.h"
#include "include/utils.h"

namespace {

constexpr unsigned int kD2Offset = 0;
constexpr unsigned int kD3Offset = 8;
constexpr unsigned int kD4Offset = 16;
constexpr int kParticleCount = 4;
constexpr int kCalibrationLayers = 5;

struct CalibrationParameters {
	double p0[kCalibrationLayers] = {0.0, 0.0, 0.0, 0.0, 0.0};
	double p1[kCalibrationLayers] = {1.0, 1.0, 1.0, 1.0, 1.0};
};

struct CutInfo {
	std::string particle;
	int charge = 0;
	int mass = 0;
	bool stop = false;
	std::unique_ptr<TCutG> cut;
};

struct Slice {
	int charge = 0;
	int mass = 0;
	bool stop = false;
	int first = -1;
	int second = -1;
};

struct ParticleCandidate {
	int charge = 0;
	int mass = 0;
	int layer = 0;
	int used_slices = 0;
	int d2 = -1;
	int d3 = -1;
	int d4 = -1;
	uint32_t mask = 0;
};

struct RebuildResult {
	double layer_energy[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double kinetic = 0.0;
	double total_energy = 0.0;
	ROOT::Math::XYZVector direction;
};

using CalculatorKey = std::tuple<int, int>;
using CalculatorMap = std::map<CalculatorKey, std::unique_ptr<brill::LostEnergyCalculator>>;

bool ParseParticleName(const std::string &particle, int &charge, int &mass);
int ElementCharge(const std::string &element);

int CutFileExists(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	std::ifstream fin(brill::CutFilePath(workspace, slice, particle, tail));
	return fin.good();
}

std::string ParticleElement(const std::string &particle) {
	size_t index = 0;
	while (
		index < particle.size()
		&& std::isdigit(static_cast<unsigned char>(particle[index]))
	) {
		++index;
	}
	return particle.substr(index);
}

bool ParseParticleName(const std::string &particle, int &charge, int &mass) {
	size_t index = 0;
	while (
		index < particle.size()
		&& std::isdigit(static_cast<unsigned char>(particle[index]))
	) {
		++index;
	}
	if (index == 0 || index >= particle.size()) return false;
	mass = std::stoi(particle.substr(0, index));
	charge = ElementCharge(particle.substr(index));
	return charge > 0;
}

int ElementCharge(const std::string &element) {
	if (element == "H") return 1;
	if (element == "He") return 2;
	if (element == "Li") return 3;
	if (element == "Be") return 4;
	if (element == "B") return 5;
	if (element == "C") return 6;
	if (element == "N") return 7;
	if (element == "O") return 8;
	return 0;
}

int BuildCut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail,
	bool required,
	CutInfo &cut
) {
	cut = CutInfo();
	cut.particle = particle;
	cut.stop = !tail;
	if (!ParseParticleName(particle, cut.charge, cut.mass)) {
		std::cerr << "Error: Parse particle name " << particle << " failed.\n";
		return -1;
	}

	std::string name = particle;
	if (!tail) {
		if (!CutFileExists(workspace, slice, particle, false)) {
			if (!required) return 1;
			std::cerr << "Error: Cut file for " << slice << "_" << particle
				<< " not exist.\n";
			return -2;
		}
	} else if (!CutFileExists(workspace, slice, particle, true)) {
		name = ParticleElement(particle);
		if (!CutFileExists(workspace, slice, name, true)) {
			if (!required) return 1;
			std::cerr << "Error: Tail cut file for " << slice << "_" << particle
				<< " not exist.\n";
			return -2;
		}
	}

	return brill::ParseCutFile(workspace, slice, name, tail, cut.cut);
}

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

void BuildSlices(
	const brill::DssdMatchEvent &left,
	const brill::DssdMatchEvent &right,
	const std::vector<CutInfo> &cuts,
	std::vector<Slice> &slices
) {
	slices.clear();
	for (int i = 0; i < left.num; ++i) {
		for (int j = 0; j < right.num; ++j) {
			for (const auto &cut : cuts) {
				if (!cut.cut || !cut.cut->IsInside(right.energy[j], left.energy[i])) continue;
				slices.push_back(Slice{
					cut.charge,
					cut.mass,
					cut.stop,
					i,
					j
				});
			}
		}
	}
}

uint32_t HitMask(int d2, int d3, int d4) {
	uint32_t mask = 0;
	if (d2 >= 0) mask |= 1u << (kD2Offset + unsigned(d2));
	if (d3 >= 0) mask |= 1u << (kD3Offset + unsigned(d3));
	if (d4 >= 0) mask |= 1u << (kD4Offset + unsigned(d4));
	return mask;
}

void BuildAlphaCandidates(
	const std::vector<Slice> &d2d3_slices,
	const std::vector<Slice> &d3d4_slices,
	std::vector<ParticleCandidate> &particles
) {
	for (const auto &slice : d2d3_slices) {
		if (slice.charge != 2 || slice.mass != 4) continue;
		if (slice.stop) {
			particles.push_back(ParticleCandidate{
				2,
				4,
				3,
				1,
				slice.first,
				slice.second,
				-1,
				HitMask(slice.first, slice.second, -1)
			});
			continue;
		}

		for (const auto &next : d3d4_slices) {
			if (next.charge != 2 || next.mass != 4) continue;
			if (next.first != slice.second) continue;
			particles.push_back(ParticleCandidate{
				2,
				4,
				next.stop ? 4 : 5,
				2,
				slice.first,
				slice.second,
				next.second,
				HitMask(slice.first, slice.second, next.second)
			});
		}
	}
}

void BuildProtonCandidates(
	const std::vector<Slice> &d2d3_slices,
	const std::vector<Slice> &d3d4_slices,
	std::vector<ParticleCandidate> &particles
) {
	for (const auto &slice : d2d3_slices) {
		if (slice.charge != 1 || slice.mass != 1) continue;
		if (slice.stop) {
			particles.push_back(ParticleCandidate{
				1,
				1,
				3,
				1,
				slice.first,
				slice.second,
				-1,
				HitMask(slice.first, slice.second, -1)
			});
			continue;
		}

		for (const auto &next : d3d4_slices) {
			if (next.charge != 1 || next.mass != 1) continue;
			if (next.first != slice.second) continue;
			particles.push_back(ParticleCandidate{
				1,
				1,
				next.stop ? 4 : 5,
				2,
				slice.first,
				slice.second,
				next.second,
				HitMask(slice.first, slice.second, next.second)
			});
		}
	}

	for (const auto &slice : d3d4_slices) {
		if (slice.charge != 1 || slice.mass != 1) continue;
		particles.push_back(ParticleCandidate{
			1,
			1,
			slice.stop ? 4 : 5,
			1,
			-1,
			slice.first,
			slice.second,
			HitMask(-1, slice.first, slice.second)
		});
	}
}

bool NonOverlapping(
	const ParticleCandidate &first,
	const ParticleCandidate &second
) {
	return (first.mask & second.mask) == 0;
}

bool SelectQuartet(
	const std::vector<ParticleCandidate> &alphas,
	const std::vector<ParticleCandidate> &protons,
	std::vector<ParticleCandidate> &selected
) {
	selected.clear();
	for (size_t i = 0; i < alphas.size(); ++i) {
		for (size_t j = i + 1; j < alphas.size(); ++j) {
			if (!NonOverlapping(alphas[i], alphas[j])) continue;
			for (size_t k = 0; k < protons.size(); ++k) {
				if (!NonOverlapping(alphas[i], protons[k])) continue;
				if (!NonOverlapping(alphas[j], protons[k])) continue;
				for (size_t l = k + 1; l < protons.size(); ++l) {
					if (!NonOverlapping(protons[k], protons[l])) continue;
					if (!NonOverlapping(alphas[i], protons[l])) continue;
					if (!NonOverlapping(alphas[j], protons[l])) continue;
					selected.push_back(alphas[i]);
					selected.push_back(alphas[j]);
					selected.push_back(protons[k]);
					selected.push_back(protons[l]);
					return true;
				}
			}
		}
	}
	return false;
}

double RecoverIncidentEnergyFromDepositedLoss(
	const brill::LostEnergyCalculator &calculator,
	double deposited_energy
) {
	if (deposited_energy <= 0.0) return 0.0;
	double maximum_loss = calculator.EnergyLossFromResidual(0.0);
	if (deposited_energy >= maximum_loss) {
		return calculator.IncidentEnergy(0.0);
	}

	double low = 0.0;
	double high = 1.0;
	while (calculator.EnergyLossFromResidual(high) > deposited_energy && high < 1e6) {
		high *= 2.0;
	}
	for (int i = 0; i < 100; ++i) {
		double middle = 0.5 * (low + high);
		double loss = calculator.EnergyLossFromResidual(middle);
		if (loss > deposited_energy) {
			low = middle;
		} else {
			high = middle;
		}
	}
	return calculator.IncidentEnergy(0.5 * (low + high));
}

const brill::LostEnergyCalculator *GetLayerCalculator(
	const brill::AppConfig &config,
	int charge,
	int mass,
	double thickness_um,
	CalculatorMap &calculators
) {
	const CalculatorKey key(charge, mass);
	auto iterator = calculators.find(key);
	if (iterator != calculators.end()) return iterator->second.get();

	auto inserted = calculators.emplace(
		key,
		std::make_unique<brill::LostEnergyCalculator>(
			charge,
			mass,
			brill::SiliconMaterial(),
			thickness_um,
			TString::Format(
				"%s/si_%.0fum_z%d_a%d.root",
				brill::JoinPath(config.root.workspace, config.paths.energy_calculator).c_str(),
				thickness_um,
				charge,
				mass
			).Data()
		)
	);
	return inserted.first->second.get();
}

bool BuildDirection(
	const ParticleCandidate &particle,
	const brill::DssdMatchEvent &d2_event,
	const brill::DssdMatchEvent &d3_event,
	ROOT::Math::XYZVector &direction
) {
	if (particle.d2 >= 0) {
		direction = ROOT::Math::XYZVector(
			d2_event.x[particle.d2],
			d2_event.y[particle.d2],
			d2_event.z[particle.d2]
		);
	} else if (particle.d3 >= 0) {
		direction = ROOT::Math::XYZVector(
			d3_event.x[particle.d3],
			d3_event.y[particle.d3],
			d3_event.z[particle.d3]
		);
	} else {
		return false;
	}
	if (direction.R() < 1e-6) return false;
	direction = direction.Unit();
	return true;
}

bool RebuildParticle(
	const ParticleCandidate &particle,
	const brill::DssdMatchEvent &d2_event,
	const brill::DssdMatchEvent &d3_event,
	const brill::DssdMatchEvent &d4_event,
	const CalibrationParameters &calibration,
	const brill::AppConfig &config,
	double d1_thickness_um,
	double d2_thickness_um,
	double d4_thickness_um,
	CalculatorMap &d1_calculators,
	CalculatorMap &d2_calculators,
	CalculatorMap &d4_calculators,
	RebuildResult &result
) {
	if (!BuildDirection(particle, d2_event, d3_event, result.direction)) return false;

	for (int i = 0; i < 6; ++i) result.layer_energy[i] = 0.0;
	double energy_before_d1 = 0.0;
	const double d2_energy = particle.d2 >= 0
		? CalibrateEnergy(calibration, 1, d2_event.energy[particle.d2])
		: 0.0;
	const double d3_energy = particle.d3 >= 0
		? CalibrateEnergy(calibration, 2, d3_event.energy[particle.d3])
		: 0.0;
	const double d4_energy = particle.d4 >= 0
		? CalibrateEnergy(calibration, 3, d4_event.energy[particle.d4])
		: 0.0;

	if (particle.layer == 3) {
		if (particle.d3 < 0) return false;
		result.layer_energy[2] = d3_energy;
		double energy_before_d2 = d3_energy;
		if (particle.d2 >= 0) {
			energy_before_d2 += d2_energy;
		} else {
			const auto *d2_calculator = GetLayerCalculator(
				config,
				particle.charge,
				particle.mass,
				d2_thickness_um,
				d2_calculators
			);
			if (!d2_calculator) return false;
			energy_before_d2 = d2_calculator->IncidentEnergy(energy_before_d2);
		}
		const auto *d1_calculator = GetLayerCalculator(
			config,
			particle.charge,
			particle.mass,
			d1_thickness_um,
			d1_calculators
		);
		if (!d1_calculator) return false;
		energy_before_d1 = d1_calculator->IncidentEnergy(energy_before_d2);
	} else if (particle.layer == 4) {
		if (particle.d4 < 0) return false;
		result.layer_energy[3] = d4_energy;
		double energy_before_d3 = d4_energy;
		if (particle.d3 >= 0) energy_before_d3 += d3_energy;
		double energy_before_d2 = energy_before_d3;
		if (particle.d2 >= 0) {
			energy_before_d2 += d2_energy;
		} else {
			const auto *d2_calculator = GetLayerCalculator(
				config,
				particle.charge,
				particle.mass,
				d2_thickness_um,
				d2_calculators
			);
			if (!d2_calculator) return false;
			energy_before_d2 = d2_calculator->IncidentEnergy(energy_before_d2);
		}
		const auto *d1_calculator = GetLayerCalculator(
			config,
			particle.charge,
			particle.mass,
			d1_thickness_um,
			d1_calculators
		);
		if (!d1_calculator) return false;
		energy_before_d1 = d1_calculator->IncidentEnergy(energy_before_d2);
	} else if (particle.layer == 5) {
		if (particle.d4 < 0) return false;
		result.layer_energy[3] = d4_energy;
		const auto *d4_calculator = GetLayerCalculator(
			config,
			particle.charge,
			particle.mass,
			d4_thickness_um,
			d4_calculators
		);
		if (!d4_calculator) return false;
		double energy_before_d3 = RecoverIncidentEnergyFromDepositedLoss(
			*d4_calculator,
			d4_energy
		);
		if (particle.d3 >= 0) energy_before_d3 += d3_energy;
		double energy_before_d2 = energy_before_d3;
		if (particle.d2 >= 0) {
			energy_before_d2 += d2_energy;
		} else {
			const auto *d2_calculator = GetLayerCalculator(
				config,
				particle.charge,
				particle.mass,
				d2_thickness_um,
				d2_calculators
			);
			if (!d2_calculator) return false;
			energy_before_d2 = d2_calculator->IncidentEnergy(energy_before_d2);
		}
		const auto *d1_calculator = GetLayerCalculator(
			config,
			particle.charge,
			particle.mass,
			d1_thickness_um,
			d1_calculators
		);
		if (!d1_calculator) return false;
		energy_before_d1 = d1_calculator->IncidentEnergy(energy_before_d2);
	} else {
		return false;
	}

	double rest_mass = brill::GetMass(particle.charge, particle.mass);
	result.kinetic = energy_before_d1;
	result.total_energy = rest_mass + energy_before_d1;
	return energy_before_d1 > 0.0;
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options(
		"rebuild_2alpha_2p",
		"Rebuild fixed 2alpha-2proton events from T0 match events."
	);
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End of run number.", cxxopts::value<int>(), "run")
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
	brill::SetAssetsPath(config.root.assets);
	const std::string calibration_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.root.workspace, config.paths.calibration).c_str()
	).Data();
	CalibrationParameters calibration;
	if (ReadCalibrationParameters(calibration_path, calibration)) {
		return 1;
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const auto *d1_detector = config.FindDetector("t0d1");
	if (!d1_detector) {
		std::cerr << "Error: Missing detector config for t0d1.\n";
		return 1;
	}
	const auto *d2_detector = config.FindDetector("t0d2");
	if (!d2_detector) {
		std::cerr << "Error: Missing detector config for t0d2.\n";
		return 1;
	}
	const auto *d4_detector = config.FindDetector("t0d4");
	if (!d4_detector) {
		std::cerr << "Error: Missing detector config for t0d4.\n";
		return 1;
	}

	std::vector<CutInfo> d2d3_cuts;
	std::vector<CutInfo> d3d4_cuts;
	CutInfo cut;
	if (BuildCut(config.root.workspace, "t0d2d3", "4He", false, true, cut)) return 1;
	d2d3_cuts.push_back(std::move(cut));
	if (BuildCut(config.root.workspace, "t0d2d3", "4He", true, true, cut)) return 1;
	d2d3_cuts.push_back(std::move(cut));
	int proton_d2d3_stop = BuildCut(config.root.workspace, "t0d2d3", "1H", false, false, cut);
	if (proton_d2d3_stop < 0) return 1;
	if (proton_d2d3_stop == 0) d2d3_cuts.push_back(std::move(cut));
	int proton_d2d3_tail = BuildCut(config.root.workspace, "t0d2d3", "1H", true, false, cut);
	if (proton_d2d3_tail < 0) return 1;
	if (proton_d2d3_tail == 0) d2d3_cuts.push_back(std::move(cut));

	if (BuildCut(config.root.workspace, "t0d3d4", "4He", false, true, cut)) return 1;
	d3d4_cuts.push_back(std::move(cut));
	if (BuildCut(config.root.workspace, "t0d3d4", "4He", true, true, cut)) return 1;
	d3d4_cuts.push_back(std::move(cut));
	if (BuildCut(config.root.workspace, "t0d3d4", "1H", false, true, cut)) return 1;
	d3d4_cuts.push_back(std::move(cut));
	if (BuildCut(config.root.workspace, "t0d3d4", "1H", true, true, cut)) return 1;
	d3d4_cuts.push_back(std::move(cut));

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);

	TChain chain2("tree");
	TChain chain3("tree");
	TChain chain4("tree");
	std::vector<int> source_runs;
	for (int current = run; current <= end_run; ++current) {
		if (config.IsSkipRun(current)) continue;
		std::string d2_path = TString::Format(
			"%s/t0d2_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current
		).Data();
		std::string d3_path = TString::Format(
			"%s/t0d3_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current
		).Data();
		std::string d4_path = TString::Format(
			"%s/t0d4_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current
		).Data();
		if (
			!std::filesystem::exists(d2_path)
			|| !std::filesystem::exists(d3_path)
			|| !std::filesystem::exists(d4_path)
		) {
			std::cerr << "Warning: Skip missing match files for run " << current << ".\n";
			continue;
		}
		chain2.Add(d2_path.c_str());
		chain3.Add(d3_path.c_str());
		chain4.Add(d4_path.c_str());
		source_runs.push_back(current);
	}
	if (source_runs.empty()) {
		std::cout << "No matched T0 files to process.\n";
		return 0;
	}
	chain2.AddFriend(&chain3, "d3");
	chain2.AddFriend(&chain4, "d4");

	brill::DssdMatchEvent d2_event;
	brill::DssdMatchEvent d3_event;
	brill::DssdMatchEvent d4_event;
	brill::SetupInput(&chain2, d2_event);
	brill::SetupInput(&chain2, d3_event, "d3.");
	brill::SetupInput(&chain2, d4_event, "d4.");

	std::string output_path = TString::Format(
		"%s/t0_2a2p_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.particle).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	).Data();
	std::filesystem::path output_file_path(output_path);
	if (!output_file_path.parent_path().empty()) {
		std::filesystem::create_directories(output_file_path.parent_path());
	}
	TFile opf(output_path.c_str(), "recreate");
	TTree opt("tree", "rebuilt 2alpha 2proton particles");
	int charge[kParticleCount] = {0, 0, 0, 0};
	int mass[kParticleCount] = {0, 0, 0, 0};
	int layer[kParticleCount] = {0, 0, 0, 0};
	double layer_energy[kParticleCount][6] = {0.0, 0.0, 0.0, 0.0};
	double kinetic[kParticleCount] = {0.0, 0.0, 0.0, 0.0};
	double energy[kParticleCount] = {0.0, 0.0, 0.0, 0.0};
	double px[kParticleCount] = {0.0, 0.0, 0.0, 0.0};
	double py[kParticleCount] = {0.0, 0.0, 0.0, 0.0};
	double pz[kParticleCount] = {0.0, 0.0, 0.0, 0.0};
	double excitation = 0.0;
	int orig_run = 0;
	long long orig_entry = -1;
	opt.Branch("charge", charge, "Z[4]/I");
	opt.Branch("mass", mass, "A[4]/I");
	opt.Branch("layer", layer, "layer[4]/I");
	opt.Branch("layer_energy", layer_energy, "le[4][6]/D");
	opt.Branch("kinetic", kinetic, "k[4]/D");
	opt.Branch("energy", energy, "e[4]/D");
	opt.Branch("px", px, "px[4]/D");
	opt.Branch("py", py, "py[4]/D");
	opt.Branch("pz", pz, "pz[4]/D");
	opt.Branch("ex", &excitation, "ex/D");
	opt.Branch("run", &orig_run, "run/I");
	opt.Branch("entry", &orig_entry, "entry/L");

	CalculatorMap d1_calculators;
	CalculatorMap d2_calculators;
	CalculatorMap d4_calculators;
	std::vector<Slice> d2d3_slices;
	std::vector<Slice> d3d4_slices;
	std::vector<ParticleCandidate> alphas;
	std::vector<ParticleCandidate> protons;
	std::vector<ParticleCandidate> quartet;

	const long long total = chain2.GetEntries();
	long long last_percentage = -1;
	std::printf("Rebuilding 2alpha2p   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}

		std::fill_n(charge, kParticleCount, 0);
		std::fill_n(mass, kParticleCount, 0);
		std::fill_n(layer, kParticleCount, 0);
		for (int i = 0; i < kParticleCount; ++i) {
			for (int j = 0; j < 6; ++j) {
				layer_energy[i][j] = 0.0;
			}
		}
		std::fill_n(kinetic, kParticleCount, 0.0);
		std::fill_n(energy, kParticleCount, 0.0);
		std::fill_n(px, kParticleCount, 0.0);
		std::fill_n(py, kParticleCount, 0.0);
		std::fill_n(pz, kParticleCount, 0.0);
		excitation = 0.0;

		long long local_entry = chain2.LoadTree(entry);
		if (local_entry < 0) continue;
		chain2.GetEntry(entry);

		BuildSlices(d2_event, d3_event, d2d3_cuts, d2d3_slices);
		BuildSlices(d3_event, d4_event, d3d4_cuts, d3d4_slices);

		alphas.clear();
		protons.clear();
		BuildAlphaCandidates(d2d3_slices, d3d4_slices, alphas);
		BuildProtonCandidates(d2d3_slices, d3d4_slices, protons);
		if (!SelectQuartet(alphas, protons, quartet)) continue;

		bool valid = true;
		for (int i = 0; i < kParticleCount; ++i) {
			const auto &particle = quartet[i];
			RebuildResult rebuild;
			if (!RebuildParticle(
				particle,
				d2_event,
				d3_event,
				d4_event,
				calibration,
				config,
				d1_detector->thickness_um,
				d2_detector->thickness_um,
				d4_detector->thickness_um,
				d1_calculators,
				d2_calculators,
				d4_calculators,
				rebuild
			)) {
				valid = false;
				break;
			}
			charge[i] = particle.charge;
			mass[i] = particle.mass;
			layer[i] = particle.layer;
			for (int j = 0; j < 5; ++j) {
				layer_energy[i][j] = rebuild.layer_energy[j];
			}
			kinetic[i] = rebuild.kinetic;
			energy[i] = rebuild.total_energy;
			px[i] = rebuild.direction.X();
			py[i] = rebuild.direction.Y();
			pz[i] = rebuild.direction.Z();
		}
		if (!valid) continue;

		brill::Particle fragments[kParticleCount] = {
			brill::Particle(
				charge[0],
				mass[0],
				kinetic[0],
				ROOT::Math::XYZVector(px[0], py[0], pz[0])
			),
			brill::Particle(
				charge[1],
				mass[1],
				kinetic[1],
				ROOT::Math::XYZVector(px[1], py[1], pz[1])
			),
			brill::Particle(
				charge[2],
				mass[2],
				kinetic[2],
				ROOT::Math::XYZVector(px[2], py[2], pz[2])
			),
			brill::Particle(
				charge[3],
				mass[3],
				kinetic[3],
				ROOT::Math::XYZVector(px[3], py[3], pz[3])
			)
		};
		brill::Particle carbon10 = fragments[0] + fragments[1];
		carbon10 = carbon10 + fragments[2];
		carbon10 = carbon10 + fragments[3];
		excitation = carbon10.ExcitationEnergy();

		int tree_number = chain2.GetTreeNumber();
		if (tree_number < 0 || size_t(tree_number) >= source_runs.size()) continue;
		orig_run = source_runs[tree_number];
		orig_entry = entry - chain2.GetTreeOffset()[tree_number];
		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	return 0;
}
