#include "include/t0/dssd.h"

#include <algorithm>
#include <cmath>

namespace brill {

namespace {

struct MatchCandidate {
	int front_index = -1;
	int back_index = -1;
	double front_strip = 0.0;
	double back_strip = 0.0;
	double energy = 0.0;
	double time = 0.0;
	int merge_tag = 0;
};

double StripPosition(double center, double size, int strips, double strip) {
	return center + size * ((strip + 0.5) / double(strips) - 0.5);
}

int SearchAdjacentFront(
	const DssdEvent &event,
	int seed_index,
	const bool *front_used
) {
	int found = -1;
	for (int i = 0; i < event.front_num; ++i) {
		if (i == seed_index) continue;
		if (front_used[i]) continue;
		if (std::abs(event.front_strip[i] - event.front_strip[seed_index]) != 1) continue;
		if (found == -1 || event.front_energy[i] > event.front_energy[found]) {
			found = i;
		}
	}
	return found;
}

int SearchAdjacentBack(
	const DssdEvent &event,
	int seed_index,
	const bool *back_used
) {
	int found = -1;
	for (int i = 0; i < event.back_num; ++i) {
		if (i == seed_index) continue;
		if (back_used[i]) continue;
		if (std::abs(event.back_strip[i] - event.back_strip[seed_index]) != 1) continue;
		if (found == -1 || event.back_energy[i] > event.back_energy[found]) {
			found = i;
		}
	}
	return found;
}

double WeightedStrip(int strip0, double energy0, int strip1, double energy1) {
	if (energy0 * 0.1 > energy1) return double(strip0);
	return 0.5 * double(strip0 + strip1);
}

void FillPhysicalPosition(
	const SiliconDetectorConfig &detector,
	double front_strip,
	double back_strip,
	double &x,
	double &y,
	double &z
) {
	if (detector.name == "t0d2") {
		x = StripPosition(
			detector.center_x_mm,
			detector.size_x_mm,
			detector.back_strips,
			back_strip
		);
		y = StripPosition(
			detector.center_y_mm,
			detector.size_y_mm,
			detector.front_strips,
			63-front_strip
		);
	} else {
		x = StripPosition(
			detector.center_x_mm,
			detector.size_x_mm,
			detector.front_strips,
			front_strip
		);
		y = StripPosition(
			detector.center_y_mm,
			detector.size_y_mm,
			detector.back_strips,
			back_strip
		);
	}
	z = detector.z_mm;
}

void AppendMatch(
	DssdMatchEvent &output,
	const SiliconDetectorConfig &detector,
	const MatchCandidate &candidate
) {
	if (output.num >= 8) return;
	int index = output.num++;
	output.front_strip[index] = candidate.front_index;
	output.back_strip[index] = candidate.back_index;
	output.energy[index] = candidate.energy;
	output.time[index] = candidate.time;
	output.merge_tag[index] = candidate.merge_tag;
	FillPhysicalPosition(
		detector,
		candidate.front_strip,
		candidate.back_strip,
		output.x[index],
		output.y[index],
		output.z[index]
	);
}

} // namespace

DssdNormalizeParameters::DssdNormalizeParameters(
	const int front_strips,
	const int back_strips
) {
	front_strips_ = front_strips;
	back_strips_ = back_strips;
}

int DssdNormalizeParameters::Write(const std::string &path) const {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file " << path << " failed.\n";
		return -1;
	}
	fout << "side strip p0 p1 p2\n";
	// front
	for (int i = 0; i < front_strips_; ++i) {
		fout << "0 " << i << " "
			<< front_p0_[i] << " "
			<< front_p1_[i] << " "
			<< front_p2_[i] << "\n";
	}
	// back
	for (int i = 0; i < back_strips_; ++i) {
		fout << "1 " << i << " "
			<< back_p0_[i] << " "
			<< back_p1_[i] << " "
			<< back_p2_[i] << "\n";
	}
	return 0;
}

int DssdNormalizeParameters::Read(const std::string &path) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file " << path << " failed.\n";
		return -1;
	}

	std::string line;
	if (!std::getline(fin, line)) {
		std::cerr << "Error: Read header from " << path << " failed.\n";
		return -1;
	}

	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		std::istringstream iss(line);
		int side = -1;
		int strip = -1;
		double value0 = 0.0;
		double value1 = 0.0;
		double value2 = 0.0;
		if (!(iss >> side >> strip >> value0 >> value1 >> value2)) continue;
		if (strip < 0 || strip >= (side == 0 ? front_strips_ : back_strips_)) {
			std::cerr << "Error: Side " << side << " strip " << strip
				<< " out of range in " << path << ".\n";
			return -1;
		}
		if (side == 0) {
			front_p0_[strip] = value0;
			front_p1_[strip] = value1;
			front_p2_[strip] = value2;
		} else {
			back_p0_[strip] = value0;
			back_p1_[strip] = value1;
			back_p2_[strip] = value2;
		}
	}
	return 0;
}


struct Hit {
	int strip = 0;
	double energy = 0.0;
	double time = 0.0;
};

void SortHits(const int num, int *strip, double *energy, double *time) {
	if (num <= 0 || num > 8) return;
	Hit hits[16];
	for (int i = 0; i < num; ++i) {
		hits[i].strip = strip[i];
		hits[i].energy = energy[i];
		hits[i].time = time[i];
	}
	std::sort(
		hits,
		hits+num,
		[](const Hit &left, const Hit &right) {
			return left.energy > right.energy;
		}
	);
	for (int i = 0; i < num; ++i) {
		strip[i] = hits[i].strip;
		energy[i] = hits[i].energy;
		time[i] = hits[i].time;
	}
}

void DssdNormalizeParameters::Apply(
	const DssdEvent &input,
	DssdEvent &output
) const {
	output.front_num = input.front_num;
	output.back_num = input.back_num;
	for (int i = 0; i < input.front_num; ++i) {
		int strip = input.front_strip[i];
		output.front_strip[i] = strip;
		output.front_energy[i] = NormEnergy(
			0, strip, input.front_energy[i]
		);
		output.front_time[i] = input.front_time[i];
	}
	for (int i = 0; i < input.back_num; ++i) {
		int strip = input.back_strip[i];
		output.back_strip[i] = strip;
		output.back_energy[i] = NormEnergy(
			1, strip, input.back_energy[i]
		);
		output.back_time[i] = input.back_time[i];
	}
	SortHits(output.front_num, output.front_strip, output.front_energy, output.front_time);
	SortHits(output.back_num, output.back_strip, output.back_energy, output.back_time);
}

void MatchDssdEvent(
	const DssdEvent &input,
	const SiliconDetectorConfig &detector,
	DssdMatchEvent &output
) {
	Reset(output);
	bool front_used[8] = {false};
	bool back_used[8] = {false};

	for (int i = 0; i < input.front_num && output.num < 8; ++i) {
		if (front_used[i]) continue;
		int fi = SearchAdjacentFront(input, i, front_used);
		bool matched = false;

		if (fi >= 0 && output.num + 1 < 8) {
			double front_total = input.front_energy[i] + input.front_energy[fi];
			for (int j = 0; j < input.back_num && !matched; ++j) {
				if (back_used[j]) continue;
				int bi = SearchAdjacentBack(input, j, back_used);
				if (bi < 0) continue;
				double back_total = input.back_energy[j] + input.back_energy[bi];
				if (std::fabs(front_total - back_total) >= detector.match_tolerance) continue;

				double diff00 = std::fabs(input.front_energy[i] - input.back_energy[j]);
				double diff11 = std::fabs(input.front_energy[fi] - input.back_energy[bi]);
				double diff01 = std::fabs(input.front_energy[i] - input.back_energy[bi]);
				double diff10 = std::fabs(input.front_energy[fi] - input.back_energy[j]);
				bool straight = diff00 < detector.match_tolerance && diff11 < detector.match_tolerance;
				bool cross = diff01 < detector.match_tolerance && diff10 < detector.match_tolerance;
				if (straight || cross) {
					bool use_cross = false;
					if (straight && cross) {
						use_cross = diff01 + diff10 < diff00 + diff11;
					} else if (cross) {
						use_cross = true;
					}
					int first_back = use_cross ? bi : j;
					int second_back = use_cross ? j : bi;

					AppendMatch(
						output,
						detector,
						MatchCandidate{
							input.front_strip[i],
							input.back_strip[first_back],
							double(input.front_strip[i]),
							double(input.back_strip[first_back]),
							input.front_energy[i],
							input.front_time[i],
							0
						}
					);
					AppendMatch(
						output,
						detector,
						MatchCandidate{
							input.front_strip[fi],
							input.back_strip[second_back],
							double(input.front_strip[fi]),
							double(input.back_strip[second_back]),
							input.front_energy[fi],
							input.front_time[fi],
							0
						}
					);
					front_used[i] = true;
					front_used[fi] = true;
					back_used[j] = true;
					back_used[bi] = true;
					matched = true;
					break;
				}

				int leading_back = input.back_energy[j] >= input.back_energy[bi] ? j : bi;
				AppendMatch(
					output,
					detector,
					MatchCandidate{
						input.front_strip[i],
						input.back_strip[leading_back],
						WeightedStrip(
							input.front_strip[i],
							input.front_energy[i],
							input.front_strip[fi],
							input.front_energy[fi]
						),
						WeightedStrip(
							input.back_strip[j],
							input.back_energy[j],
							input.back_strip[bi],
							input.back_energy[bi]
						),
						front_total,
						input.front_time[i],
						3
					}
				);
				front_used[i] = true;
				front_used[fi] = true;
				back_used[j] = true;
				back_used[bi] = true;
				matched = true;
			}
		}

		if (!matched && fi >= 0) {
			double front_total = input.front_energy[i] + input.front_energy[fi];
			for (int j = 0; j < input.back_num && !matched; ++j) {
				if (back_used[j]) continue;
				if (std::fabs(front_total - input.back_energy[j]) >= detector.match_tolerance) continue;
				AppendMatch(
					output,
					detector,
					MatchCandidate{
						input.front_strip[i],
						input.back_strip[j],
						WeightedStrip(
							input.front_strip[i],
							input.front_energy[i],
							input.front_strip[fi],
							input.front_energy[fi]
						),
						double(input.back_strip[j]),
						front_total,
						input.front_time[i],
						2
					}
				);
				front_used[i] = true;
				front_used[fi] = true;
				back_used[j] = true;
				matched = true;
			}
		}

		if (!matched) {
			for (int j = 0; j < input.back_num && !matched; ++j) {
				if (back_used[j]) continue;
				int bi = SearchAdjacentBack(input, j, back_used);
				if (bi < 0) continue;
				double back_total = input.back_energy[j] + input.back_energy[bi];
				if (std::fabs(input.front_energy[i] - back_total) >= detector.match_tolerance) continue;
				int leading_back = input.back_energy[j] >= input.back_energy[bi] ? j : bi;
				AppendMatch(
					output,
					detector,
					MatchCandidate{
						input.front_strip[i],
						input.back_strip[leading_back],
						double(input.front_strip[i]),
						WeightedStrip(
							input.back_strip[j],
							input.back_energy[j],
							input.back_strip[bi],
							input.back_energy[bi]
						),
						input.front_energy[i],
						input.front_time[i],
						1
					}
				);
				front_used[i] = true;
				back_used[j] = true;
				back_used[bi] = true;
				matched = true;
			}
		}

		if (!matched) {
			for (int j = 0; j < input.back_num && !matched; ++j) {
				if (back_used[j]) continue;
				if (std::fabs(input.front_energy[i] - input.back_energy[j]) >= detector.match_tolerance) continue;
				AppendMatch(
					output,
					detector,
					MatchCandidate{
						input.front_strip[i],
						input.back_strip[j],
						double(input.front_strip[i]),
						double(input.back_strip[j]),
						input.front_energy[i],
						input.front_time[i],
						0
					}
				);
				front_used[i] = true;
				back_used[j] = true;
				matched = true;
			}
		}
	}
}

} // namespace brill
