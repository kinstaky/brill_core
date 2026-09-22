#ifndef __BRILL_DELTA_ENERGY_CALCULATOR_H__
#define __BRILL_DELTA_ENERGY_CALCULATOR_H__

#include <memory>
#include <string>
#include <vector>

#include <TSpline.h>

#include "include/config.h"

namespace brill {

class DeltaEnergyCalculator {
public:
	DeltaEnergyCalculator(
		const AppConfig &config,
		int charge,
		int mass,
		int max_energy = 1000
	);

	double Energy(unsigned short layer, const double delta_energy) const;
	double DeltaEnergy(unsigned short layer, const double energy) const;

	int Initialize(const AppConfig &config);

private:
	int charge_ = 0;
	int mass_ = 0;
	int max_energy_ = 1000;
	std::vector<double> thickness_;
	std::vector<std::unique_ptr<TSpline3>> de_e_funcs_;
	std::vector<std::unique_ptr<TSpline3>> e_de_funcs_;

	int Load(const AppConfig &config);
	std::string CachePath(const AppConfig &config) const;
};

} // namespace brill

#endif
