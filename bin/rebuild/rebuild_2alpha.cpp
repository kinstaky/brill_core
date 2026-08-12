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

bool SelectTwoAlpha(
	const brill::ParticleEvent &event,
	int &first,
	int &second
) {
	first = -1;
	second = -1;
	for (int i = 0; i < event.num; ++i) {
		if (event.charge[i] != 2 || event.mass[i] != 4) continue;
		if (first < 0) {
			first = i;
		} else if (second < 0) {
			second = i;
		} else {
			return false;
		}
	}
	return first >= 0 && second >= 0;
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
	cxxopts::Options options("rebuild_2alpha", "Rebuild 2-alpha spectrum from T0 particles.");
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
		if (config.IsJumpRun(current_run)) continue;
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
		"%s/Be8_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.spectrum).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	).Data();
	std::filesystem::path output_file_path(output_path);
	TFile opf(output_path.c_str(), "recreate");
	TH1F ex("ex", "8Be excitation energy", 600, -5.0, 25.0);
	TTree opt("tree", "rebuilt 8Be");
	double excitation;
	double energy[2];
	double px[2];
	double py[2];
	double pz[2];
	double x[2];
	double y[2];
	double z[2];
	opt.Branch("excitation", &excitation, "ex/D");
	opt.Branch("energy", energy, "e[2]/D");
	opt.Branch("px", px, "px[2]/D");
	opt.Branch("py", py, "py[2]/D");
	opt.Branch("pz", pz, "pz[2]/D");
	opt.Branch("x", x, "x[2]/D");
	opt.Branch("y", y, "y[2]/D");
	opt.Branch("z", z, "z[2]/D");


	const long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Rebuilding 2alpha   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain.GetEntry(entry);

		int first = -1;
		int second = -1;
		if (!SelectTwoAlpha(event, first, second)) continue;

		ROOT::Math::XYZVector direction0;
		ROOT::Math::XYZVector direction1;
		if (!BuildDirection(event, first, direction0)) continue;
		if (!BuildDirection(event, second, direction1)) continue;

		double energy0 = TotalAlphaEnergy(event, first, d1_calculator, downstream_calculator);
		double energy1 = TotalAlphaEnergy(event, second, d1_calculator, downstream_calculator);
		if (energy0 <= 0.0 || energy1 <= 0.0) continue;

		brill::Particle alpha0(2, 4, energy0, direction0);
		brill::Particle alpha1(2, 4, energy1, direction1);
		brill::Particle be8 = alpha0 + alpha1;
		ex.Fill(be8.ExcitationEnergy());

		excitation = be8.ExcitationEnergy();
		energy[0] = alpha0.Energy();
		energy[1] = alpha1.Energy();
		px[0] = alpha0.Direction().X();
		px[1] = alpha1.Direction().X();
		py[0] = alpha0.Direction().Y();
		py[1] = alpha1.Direction().Y();
		pz[0] = alpha0.Direction().Z();
		pz[1] = alpha1.Direction().Z();
		x[0] = event.x[first];
		x[1] = event.x[second];
		y[0] = event.y[first];
		y[1] = event.y[second];
		z[0] = event.z[first];
		z[1] = event.z[second];
		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	ex.Write();
	opt.Write();
	opf.Close();
	return 0;
}
