#include "include/energy_calculator/delta_energy_calculator.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TString.h>

#include "include/energy_calculator/range_energy_calculator.h"
#include "include/utils.h"

namespace {

constexpr double kEnergyStep = 0.1;

bool ComparePairFirst(
	const std::pair<double, double> &left,
	const std::pair<double, double> &right
) {
	return left.first < right.first;
}

int CollectSiliconDetectors(
	const brill::AppConfig &config,
	std::vector<double> &thickness
) {
	thickness.clear();
	for (const auto &name : config.t0_.silicon) {
		const auto *detector = config.FindDetector(name);
		if (!detector) {
			std::cerr << "Error: T0 silicon detector " << name << " not found in config.\n";
			return -1;
		}
		thickness.push_back(detector->thickness_um);
	}
	if (thickness.size() < 2) {
		std::cerr << "Error: Need at least 2 T0 silicon detectors in config.t0.silicon.\n";
		return -1;
	}
	return 0;
}

std::string RangeCachePath(
	const brill::AppConfig &config,
	int charge,
	int mass
) {
	return TString::Format(
		"%s/si_z%d_a%d.root",
		brill::JoinPath(config.root_.workspace, config.paths_.energy_calculator).c_str(),
		charge,
		mass
	).Data();
}

int BuildSliceFunctions(
	const brill::RangeEnergyCalculator &calculator,
	double first_thickness,
	double second_thickness,
	std::vector<double> &delta_energy,
	std::vector<double> &energy
) {
	delta_energy.clear();
	energy.clear();

	double min_total_energy = calculator.Energy(first_thickness);
	double max_total_energy = calculator.Energy(first_thickness + second_thickness);

	for (
		double total_energy = min_total_energy;
		total_energy <= max_total_energy;
		total_energy += kEnergyStep
	) {
		double residual_range = calculator.Range(total_energy) - first_thickness;
		if (residual_range <= 0.0) continue;
		if (residual_range > second_thickness) break;

		double residual_energy = calculator.Energy(residual_range);
		if (residual_energy < 0.0 || residual_energy > total_energy) continue;

		delta_energy.push_back(total_energy - residual_energy);
		energy.push_back(residual_energy);
	}

	return delta_energy.size() >= 3 ? 0 : -1;
}

} // namespace

namespace brill {

DeltaEnergyCalculator::DeltaEnergyCalculator(
	const AppConfig &config,
	int charge,
	int mass
)
: charge_(charge)
, mass_(mass) {
	if (Load(config) != 0) {
		if (Initialize(config, charge_, mass_)) {
			throw std::runtime_error("Initialize T0 delta energy calculator failed.");
		}
		if (Load(config) != 0) {
			throw std::runtime_error("Load T0 delta energy calculator failed.");
		}
	}
}

double DeltaEnergyCalculator::Energy(unsigned short layer, double delta_energy) const {
	if (layer >= e_de_funcs_.size() || !e_de_funcs_[layer]) return 0.0;
	return e_de_funcs_[layer]->Eval(delta_energy);
}

double DeltaEnergyCalculator::DeltaEnergy(unsigned short layer, double energy) const {
	if (layer >= de_e_funcs_.size() || !de_e_funcs_[layer]) return 0.0;
	return de_e_funcs_[layer]->Eval(energy);
}

int DeltaEnergyCalculator::Initialize(
	const AppConfig &config,
	int charge,
	int mass
) {
	if (CollectSiliconDetectors(config, thickness_)) return -1;

	std::filesystem::path path(
		TString::Format(
			"%s/t0_delta_z%d_a%d.root",
			JoinPath(config.root_.workspace, config.paths_.energy_calculator).c_str(),
			charge,
			mass
		).Data()
	);
	if (!path.parent_path().empty()) {
		std::filesystem::create_directories(path.parent_path());
	}

	RangeEnergyCalculator calculator(
		charge,
		mass,
		SiliconMaterial(),
		RangeCachePath(config, charge, mass)
	);

	TFile output(path.string().c_str(), "recreate");
	if (output.IsZombie()) {
		std::cerr << "Error: Open " << path.string() << " failed.\n";
		return -1;
	}

	std::vector<double> delta_energy;
	std::vector<double> energy;
	for (size_t i = 0; i + 1 < thickness_.size(); ++i) {
		if (BuildSliceFunctions(
			calculator,
			thickness_[i],
			thickness_[i+1],
			delta_energy,
			energy
		)) {
			std::cerr << "Error: Build delta-energy spline for T0 layer " << i << " failed.\n";
			return -1;
		}

		std::vector<std::pair<double, double>> de_e_pairs;
		std::vector<std::pair<double, double>> e_de_pairs;
		for (size_t j = 0; j < energy.size(); ++j) {
			de_e_pairs.push_back(std::make_pair(energy[j], delta_energy[j]));
			e_de_pairs.push_back(std::make_pair(delta_energy[j], energy[j]));
		}
		std::sort(de_e_pairs.begin(), de_e_pairs.end(), ComparePairFirst);
		std::sort(e_de_pairs.begin(), e_de_pairs.end(), ComparePairFirst);

		std::vector<double> de_e_x;
		std::vector<double> de_e_y;
		std::vector<double> e_de_x;
		std::vector<double> e_de_y;
		for (const auto &point : de_e_pairs) {
			de_e_x.push_back(point.first);
			de_e_y.push_back(point.second);
		}
		for (const auto &point : e_de_pairs) {
			e_de_x.push_back(point.first);
			e_de_y.push_back(point.second);
		}

		TSpline3 de_e_func(
			TString::Format("de_e_%zu", i),
			de_e_x.data(),
			de_e_y.data(),
			int(de_e_x.size())
		);
		TSpline3 e_de_func(
			TString::Format("e_de_%zu", i),
			e_de_x.data(),
			e_de_y.data(),
			int(e_de_x.size())
		);
		output.cd();
		de_e_func.Write(TString::Format("de_e_%zu", i));
		e_de_func.Write(TString::Format("e_de_%zu", i));
	}

	return 0;
}

int DeltaEnergyCalculator::Load(const AppConfig &config) {
	if (CollectSiliconDetectors(config, thickness_)) return -1;

	de_e_funcs_.clear();
	e_de_funcs_.clear();

	std::unique_ptr<TFile> file(TFile::Open(CachePath(config).c_str(), "read"));
	if (!file || file->IsZombie()) return -1;

	for (size_t i = 0; i + 1 < thickness_.size(); ++i) {
		auto *de_e = dynamic_cast<TSpline3*>(file->Get(TString::Format("de_e_%zu", i)));
		auto *e_de = dynamic_cast<TSpline3*>(file->Get(TString::Format("e_de_%zu", i)));
		if (!de_e || !e_de) return -1;
		de_e_funcs_.push_back(std::unique_ptr<TSpline3>(
			dynamic_cast<TSpline3*>(de_e->Clone(TString::Format("de_e_%zu_clone", i)))
		));
		e_de_funcs_.push_back(std::unique_ptr<TSpline3>(
			dynamic_cast<TSpline3*>(e_de->Clone(TString::Format("e_de_%zu_clone", i)))
		));
		if (!de_e_funcs_.back() || !e_de_funcs_.back()) return -1;
	}
	return 0;
}

std::string DeltaEnergyCalculator::CachePath(const AppConfig &config) const {
	return TString::Format(
		"%s/t0_delta_z%d_a%d.root",
		JoinPath(config.root_.workspace, config.paths_.energy_calculator).c_str(),
		charge_,
		mass_
	).Data();
}

} // namespace brill
