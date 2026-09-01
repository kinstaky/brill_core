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

inline bool ParseCutVectorLines(
	const std::vector<std::string> &lines,
	size_t &index,
	std::vector<double> &values
) {
	size_t left = lines[index].find('{');
	size_t right = lines[index].rfind('}');
	if (left == std::string::npos) return false;
	// in single line
	if (right != std::string::npos) {
		std::string text = lines[index].substr(left + 1, right - left - 1);
		for (char &ch : text) {
			if (ch == ',') ch = ' ';
		}
		std::istringstream iss(text);
		double value = 0.0;
		values.clear();
		while (iss >> value) {
			values.push_back(value);
		}
		index++;
		return !values.empty();
	}
	// in multiple lines
	values.clear();
	++index;
	while (index < lines.size()) {
		size_t right = lines[index].rfind('}');
		if (right != std::string::npos) {
			++index;
			break;
		}
		std::string text = lines[index];
		for (char &ch : text) {
			if (ch == ',') ch = ' ';
		}
		std::istringstream iss(text);
		double value = 0.0;
		while (iss >> value) {
			values.push_back(value);
		}
		++index;
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
	size_t line_index = 4;
	ParseCutVectorLines(lines, line_index, x);
	ParseCutVectorLines(lines, line_index, y);
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
