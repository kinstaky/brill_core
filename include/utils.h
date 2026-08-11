#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <TCutG.h>

#include "include/config.h"

namespace brill {

inline std::string JoinPath(const std::string &left, const std::string &right) {
	if (left.empty()) return right;
	if (right.empty()) return left;
	if (left.back() == '/') return left + right;
	return left + "/" + right;
}

inline std::string TriggerInfix(const std::string &trigger) {
	return trigger == "main" ? "" : (trigger + "_");
}

inline std::string CutName(
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	return slice + "_" + particle + (tail ? "_tail" : "");
}

inline std::string CutFilePath(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	return JoinPath(
		JoinPath(workspace, "cuts"),
		CutName(slice, particle, tail) + ".cpp"
	);
}

inline bool ParseCutVectorLine(const std::string &line, std::vector<double> &values) {
	size_t left = line.find('{');
	size_t right = line.rfind('}');
	if (left == std::string::npos || right == std::string::npos || right <= left) {
		return false;
	}
	std::string text = line.substr(left + 1, right - left - 1);
	for (char &ch : text) {
		if (ch == ',') ch = ' ';
	}
	std::istringstream iss(text);
	double value = 0.0;
	values.clear();
	while (iss >> value) {
		values.push_back(value);
	}
	return !values.empty();
}

inline int ParseCutFile(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail,
	std::unique_ptr<TCutG> &cut
) {
	std::string path = CutFilePath(workspace, slice, particle, tail);
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open cut file " << path << " failed.\n";
		return -1;
	}

	std::vector<std::string> lines;
	std::string line;
	while (std::getline(fin, line)) {
		lines.push_back(line);
	}

	std::vector<double> x;
	std::vector<double> y;
	for (size_t i = 4; i < lines.size(); ++i) {
		if (x.empty()) {
			ParseCutVectorLine(lines[i], x);
			if (!x.empty()) continue;
		} else if (y.empty()) {
			ParseCutVectorLine(lines[i], y);
			if (!y.empty()) break;
		}
	}
	if (x.empty() || y.empty()) {
		std::cerr << "Error: Parse cut vectors from " << path << " failed.\n";
		return -1;
	}
	if (x.size() != y.size()) {
		std::cerr << "Error: X/Y point numbers mismatch in " << path << ".\n";
		return -1;
	}

	std::string name = CutName(slice, particle, tail);
	cut = std::make_unique<TCutG>(name.c_str(), int(x.size()), x.data(), y.data());
	cut->SetName(name.c_str());
	cut->SetTitle(name.c_str());
	return 0;
}

} // namespace brill
