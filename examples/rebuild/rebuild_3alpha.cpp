#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/lost_energy_calculator.h"
#include "include/event/particle_event.h"
#include "include/rebuild/nuclear_data.h"
#include "include/rebuild/particle.h"
#include "include/utils.h"

namespace {

constexpr int kAlphaCount = 3;

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double SumDownstreamThickness(const brill::AppConfig &config) {
	double thickness = 0.0;
	for (size_t i = 1; i < config.t0.silicon.size(); ++i) {
		const auto *detector = config.FindDetector(config.t0.silicon[i]);
		if (detector) thickness += detector->thickness_um;
	}
	return thickness;
}

double RecoverEnergyWithCsI(
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
	double residual = 0.5 * (low + high);
	return calculator.IncidentEnergy(residual);
}

void CollectAlphaIndices(
	const brill::ParticleEvent &event,
	std::vector<int> &indices
) {
	indices.clear();
	for (int i = 0; i < event.num; ++i) {
		if (event.charge[i] == 2 && event.mass[i] == 4) {
			indices.push_back(i);
		}
	}
}

bool BuildDirection(
	const brill::ParticleEvent &event,
	int index,
	ROOT::Math::XYZVector &direction
) {
	direction = ROOT::Math::XYZVector(
		event.x[index],
		event.y[index],
		event.z[index]
	);
	if (direction.R() < 1e-6) return false;
	direction = direction.Unit();
	return true;
}

double TotalAlphaEnergy(
	const brill::ParticleEvent &event,
	int index,
	const brill::LostEnergyCalculator &d1_calculator,
	const brill::LostEnergyCalculator &downstream_calculator
) {
	double energy_after_d1 = event.energy[index];
	if (!event.stop[index]) {
		energy_after_d1 = RecoverEnergyWithCsI(downstream_calculator, event.energy[index]);
	}
	return d1_calculator.IncidentEnergy(energy_after_d1);
}

} // namespace

int main(int argc, char **argv) {
	cxxopts::Options options("rebuild_3alpha", "Rebuild 3-alpha spectrum from T0 particles.");
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
	brill::SetAssetsPath(config.root.assets);

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
	double downstream_thickness = SumDownstreamThickness(config);
	if (downstream_thickness <= 0.0) {
		std::cerr << "Error: Invalid downstream T0 silicon thickness.\n";
		return 1;
	}

	brill::LostEnergyCalculator d1_calculator(
		2,
		4,
		brill::SiliconMaterial(),
		d1_detector->thickness_um,
		TString::Format(
			"%s/si_%.0fum_z2_a4.root",
			brill::JoinPath(config.root.workspace, config.paths.energy_calculator).c_str(),
			d1_detector->thickness_um
		).Data()
	);
	brill::LostEnergyCalculator downstream_calculator(
		2,
		4,
		brill::SiliconMaterial(),
		downstream_thickness,
		TString::Format(
			"%s/si_%.0fum_z2_a4.root",
			brill::JoinPath(config.root.workspace, config.paths.energy_calculator).c_str(),
			downstream_thickness
		).Data()
	);

	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);
	const std::string particle_dir = brill::JoinPath(config.root.workspace, config.paths.particle);

	TChain chain("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsSkipRun(current_run)) continue;
		std::string path = TString::Format(
			"%s/t0_%s%04d.root",
			particle_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		).Data();
		if (!std::filesystem::exists(path)) {
			std::cerr << "Warning: Skip missing rebuilt particle file " << path << ".\n";
			continue;
		}
		chain.Add(path.c_str());
		++added_runs;
	}
	if (added_runs == 0) {
		std::cout << "No rebuilt particle files to process.\n";
		return 0;
	}

	brill::ParticleEvent event;
	brill::SetupInput(&chain, event);

	std::string output_path = TString::Format(
		"%s/C12_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.spectrum).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	).Data();
	std::filesystem::path output_file_path(output_path);
	if (!output_file_path.parent_path().empty()) {
		std::filesystem::create_directories(output_file_path.parent_path());
	}
	TFile opf(output_path.c_str(), "recreate");
	TH1F ex("ex", "12C excitation energy", 200, 7.0, 37.0);
	TTree opt("tree", "rebuilt 12C");
	double excitation = 0.0;
	double energy[kAlphaCount] = {0.0, 0.0, 0.0};
	double px[kAlphaCount] = {0.0, 0.0, 0.0};
	double py[kAlphaCount] = {0.0, 0.0, 0.0};
	double pz[kAlphaCount] = {0.0, 0.0, 0.0};
	double x[kAlphaCount] = {0.0, 0.0, 0.0};
	double y[kAlphaCount] = {0.0, 0.0, 0.0};
	double z[kAlphaCount] = {0.0, 0.0, 0.0};
	opt.Branch("excitation", &excitation, "ex/D");
	opt.Branch("energy", energy, "e[3]/D");
	opt.Branch("px", px, "px[3]/D");
	opt.Branch("py", py, "py[3]/D");
	opt.Branch("pz", pz, "pz[3]/D");
	opt.Branch("x", x, "x[3]/D");
	opt.Branch("y", y, "y[3]/D");
	opt.Branch("z", z, "z[3]/D");

	std::vector<int> alpha_indices;
	const long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Rebuilding 3alpha   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain.GetEntry(entry);

		CollectAlphaIndices(event, alpha_indices);
		if (alpha_indices.size() < kAlphaCount) continue;

		for (size_t i = 0; i + 2 < alpha_indices.size(); ++i) {
			for (size_t j = i + 1; j + 1 < alpha_indices.size(); ++j) {
				for (size_t k = j + 1; k < alpha_indices.size(); ++k) {
					int indices[kAlphaCount] = {
						alpha_indices[i],
						alpha_indices[j],
						alpha_indices[k]
					};
					ROOT::Math::XYZVector directions[kAlphaCount];
					bool valid = true;
					for (int alpha = 0; alpha < kAlphaCount; ++alpha) {
						if (!BuildDirection(event, indices[alpha], directions[alpha])) {
							valid = false;
							break;
						}
					}
					if (!valid) continue;

					double total_alpha_energy[kAlphaCount] = {0.0, 0.0, 0.0};
					for (int alpha = 0; alpha < kAlphaCount; ++alpha) {
						total_alpha_energy[alpha] = TotalAlphaEnergy(
							event,
							indices[alpha],
							d1_calculator,
							downstream_calculator
						);
						if (total_alpha_energy[alpha] <= 0.0) {
							valid = false;
							break;
						}
					}
					if (!valid) continue;

					brill::Particle carbon(
						2,
						4,
						total_alpha_energy[0],
						directions[0]
					);
					for (int alpha = 0; alpha < kAlphaCount; ++alpha) {
						brill::Particle current_alpha(
							2,
							4,
							total_alpha_energy[alpha],
							directions[alpha]
						);
						if (alpha > 0) {
							carbon = carbon + current_alpha;
						}
						energy[alpha] = current_alpha.Energy();
						px[alpha] = current_alpha.Direction().X();
						py[alpha] = current_alpha.Direction().Y();
						pz[alpha] = current_alpha.Direction().Z();
						x[alpha] = event.x[indices[alpha]];
						y[alpha] = event.y[indices[alpha]];
						z[alpha] = event.z[indices[alpha]];
					}

					excitation = carbon.ExcitationEnergy();
					ex.Fill(excitation);
					opt.Fill();
				}
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	ex.Write();
	opt.Write();
	opf.Close();
	return 0;
}
