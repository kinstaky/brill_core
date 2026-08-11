#include "include/config.h"

#include <iostream>

#include "external/toml.hpp"

namespace brill {

namespace {

void LoadPath(const toml::table &table, const char *key, std::string &value) {
	if (auto node = table[key].value<std::string>()) {
		value = *node;
	}
}

template<typename T>
void LoadValue(
	const toml::table &table,
	const char *key,
	T &value
) {
	if (auto node = table[key].value<T>()) {
		value = *node;
	}
}

int LoadRoot(const toml::table &table, RootConfig &root) {
	if (auto node = table["workspace"].value<std::string>()) {
		root.workspace = *node;
	} else {
		std::cerr << "Error: Missing required key workspace.\n";
		return -1;
	}
	if (auto node = table["trigger"].value<std::string>()) {
		root.trigger = *node;
	}
	if (auto node = table["assets"].value<std::string>()) {
		root.assets = *node;
	}
	if (const auto *jump_run = table["jump_run"].as_array()) {
		for (size_t i = 0; i < jump_run->size(); ++i) {
			if (auto value = (*jump_run)[i].value<int>()) {
				root.jump_run.push_back(*value);
			}
		}
	}
	return 0;
}

void LoadPaths(const toml::table &table, AppPaths &paths) {
	LoadPath(table, "ore", paths.ore);
	LoadPath(table, "grit", paths.grit);
	LoadPath(table, "ingot", paths.ingot);
	LoadPath(table, "normalize", paths.normalize);
	LoadPath(table, "match", paths.match);
	LoadPath(table, "track", paths.track);
	LoadPath(table, "particle", paths.particle);
	LoadPath(table, "estimate", paths.estimate);
	LoadPath(table, "spectrum", paths.spectrum);
	LoadPath(table, "calibration", paths.calibration);
	LoadPath(table, "energy_calculator", paths.energy_calculator);
}

template<typename T>
void LoadNormalizeStripWindow(
	const toml::table &table,
	const char *key,
	T *fill
) {
	const auto *range = table[key].as_array();
	if (!range || range->size() != 2) return;
	std::optional<T> range_start = (*range)[0].template value<T>();
	std::optional<T> range_end = (*range)[1].template value<T>();
	if (range_start) fill[0] = *range_start;
	if (range_end) fill[1] = *range_end;
}

void LoadNormalizeStrips(
	const toml::table &table,
	NromalizeStripsConfig &strip,
	int &index
) {
	LoadValue(table, "index", index);
	LoadValue(table, "norm_side", strip.norm_side);
	LoadNormalizeStripWindow(table, "ref", strip.ref);
	LoadNormalizeStripWindow(table, "norm", strip.norm);
	LoadNormalizeStripWindow(table, "ref_energy", strip.ref_energy);
	LoadNormalizeStripWindow(table, "norm_energy", strip.norm_energy);
}

void LoadNormalizeDetector(
	const toml::table &table,
	NormalizeDetectorConfig &detector
) {
	if (const auto *strips = table["strips"].as_array()) {
		detector.strips.resize(strips->size());
		for (size_t i = 0; i < strips->size(); ++i) {
			const auto *strips_table = (*strips)[i].as_table();
			if (!strips_table) continue;
			NromalizeStripsConfig strip_config;
			int index = -1;
			LoadNormalizeStrips(*strips_table, strip_config, index);
			detector.strips[index] = strip_config;
		}
	}
}

void LoadNormalize(const toml::table &table, NormalizeConfig &normalize) {
	if (const auto *run_array = table["run"].as_array()) {
		for (const auto &run_item : *run_array) {
			const auto &run_table = *run_item.as_table();
			int start_run = 0;
			int use_run = 0;
			LoadValue(run_table, "run", start_run);
			LoadValue(run_table, "use", use_run);
			normalize.runs.push_back(std::make_pair(start_run, use_run));
		}
	}
	for (const auto &item : table) {
		if (!item.second.is_table()) continue;
		NormalizeDetectorConfig detector;
		LoadNormalizeDetector(*item.second.as_table(), detector);
		normalize.detectors[std::string(item.first.str())] = detector;
	}

	std::sort(
		normalize.runs.begin(),
		normalize.runs.end(),
		[](const std::pair<int, int>& a, const std::pair<int, int> &b) {
			return a.first < b.first;
		}
	);
}
void LoadTrackWindow(
	const toml::table &table,
	const char *key,
	TrackWindowConfig &window
) {
	const auto *node = table[key].as_array();
	if (!node || node->size() < 2) return;
	auto min = (*node)[0].value<double>();
	auto max = (*node)[1].value<double>();
	if (min) window.min = *min;
	if (max) window.max = *max;
}

void LoadTrack(const toml::table &table, TrackConfig &track) {
	LoadTrackWindow(table, "d2d1_window", track.d2d1_window);
	LoadTrackWindow(table, "d3d2_window", track.d3d2_window);
	LoadTrackWindow(table, "d4d3_window", track.d4d3_window);
}

void LoadT0(const toml::table &table, T0Config &t0) {
	if (const auto *silicon = table["silicon"].as_array()) {
		t0.silicon.clear();
		for (size_t i = 0; i < silicon->size(); ++i) {
			if (auto value = (*silicon)[i].value<std::string>()) {
				t0.silicon.push_back(*value);
			}
		}
	}
}

void LoadStraightParticle(
	const toml::table &table,
	StraightParticleConfig &particle
) {
	if (auto node = table["particle"].value<std::string>()) {
		particle.particle = *node;
	}
	LoadValue(table, "mean", particle.mean);
	LoadValue(table, "sigma", particle.sigma);
	LoadValue(table, "range", particle.range);
}

void LoadStraightSlice(
	const toml::table &table,
	StraightSliceConfig &slice
) {
	LoadValue(table, "A", slice.A);
	LoadValue(table, "B", slice.B);
	if (const auto *particles = table["particles"].as_array()) {
		for (size_t i = 0; i < particles->size(); ++i) {
			const auto *particle_table = (*particles)[i].as_table();
			if (!particle_table) continue;
			StraightParticleConfig particle;
			LoadStraightParticle(*particle_table, particle);
			if (!particle.particle.empty()) {
				slice.particles.push_back(particle);
			}
		}
	}
}

void LoadIdentify(const toml::table &table, IdentifyConfig &identify) {
	const auto *straight = table["straight"].as_table();
	if (!straight) return;
	for (const auto &item : *straight) {
		if (!item.second.is_table()) continue;
		StraightSliceConfig slice;
		LoadStraightSlice(*item.second.as_table(), slice);
		identify.straight[std::string(item.first.str())] = slice;
	}
}

void LoadPpac(const toml::table &table, PpacConfig &ppac) {
	if (const auto *node = table["z_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.z_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["x_offset_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.x_offset_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["y_offset_mm"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.y_offset_mm[i] = *value;
			}
		}
	}
	if (const auto *node = table["x_scale"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.x_scale[i] = *value;
			}
		}
	}
	if (const auto *node = table["y_scale"].as_array()) {
		for (size_t i = 0; i < node->size() && i < kMaxPpac; ++i) {
			if (auto value = (*node)[i].value<double>()) {
				ppac.y_scale[i] = *value;
			}
		}
	}
}

void LoadDetector(const toml::table &table, const std::string &name, SiliconDetectorConfig &detector) {
	detector.name = name;
	if (auto node = table["type"].value<std::string>()) {
		detector.type = *node;
	}
	LoadValue(table, "front_strips", detector.front_strips);
	LoadValue(table, "back_strips", detector.back_strips);
	LoadValue(table, "thickness_um", detector.thickness_um);
	LoadValue(table, "size_x_mm", detector.size_x_mm);
	LoadValue(table, "size_y_mm", detector.size_y_mm);
	LoadValue(table, "center_x_mm", detector.center_x_mm);
	LoadValue(table, "center_y_mm", detector.center_y_mm);
	LoadValue(table, "z_mm", detector.z_mm);
	LoadValue(table, "match_tolerance", detector.match_tolerance);
	LoadValue(table, "use_integral", detector.use_integral);
}

} // namespace

int AppConfig::Load(const std::string &path) {
	try {
		*this = AppConfig();
		toml::table table = toml::parse_file(path);
		if (LoadRoot(table, root_)) return -1;
		if (const auto *paths_tbl = table["paths"].as_table()) {
			LoadPaths(*paths_tbl, paths_);
		}
		if (const auto *t0_tbl = table["t0"].as_table()) {
			LoadT0(*t0_tbl, t0_);
		}
		if (const auto *normalize_tbl = table["normalize"].as_table()) {
			LoadNormalize(*normalize_tbl, normalize_);
		}
		if (const auto *track_tbl = table["track"].as_table()) {
			LoadTrack(*track_tbl, track_);
		}
		if (const auto *identify_tbl = table["identify"].as_table()) {
			LoadIdentify(*identify_tbl, identify_);
		}

		if (const auto *detectors_tbl = table["detectors"].as_table()) {
			if (const auto *ppac_tbl = (*detectors_tbl)["ppac"].as_table()) {
				LoadPpac(*ppac_tbl, ppac_);
			}
			for (const char *name : {"t0d1", "t0d2", "t0d3", "t0d4", "t0s"}) {
				if (const auto *detector_tbl = (*detectors_tbl)[name].as_table()) {
					SiliconDetectorConfig config_detector;
					LoadDetector(*detector_tbl, name, config_detector);
					detectors_[config_detector.name] = config_detector;
				}
			}
		}
	} catch (const toml::parse_error &err) {
		std::cerr << "Error: Parsing " << path << " failed:\n" << err << "\n";
		return -1;
	}
	return 0;
}

const SiliconDetectorConfig* AppConfig::FindDetector(
	const std::string &name
) const {
	auto iter = detectors_.find(name);
	if (iter == detectors_.end()) return nullptr;
	return &iter->second;
}

const StraightSliceConfig* AppConfig::FindStraightSlice(
	const std::string &name
) const {
	auto iter = identify_.straight.find(name);
	if (iter == identify_.straight.end()) return nullptr;
	return &iter->second;
}

} // namespace brill
