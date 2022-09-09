/*
 * TSpecies.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include "BirpTypes.h"
#include <algorithm>
#include <vector>

class TDataPoint {
private:
	TypeIndex timepoint;
	TypeCounts counts;
	TypeEffort effort;

	bool operator<(const TDataPoint &other) { return timepoint < other.timepoint; }
};

class TDataSeries {
private:
	std::vector<TDataPoint> _dataPoints;

public:
	TDataSeries() = default;
	size_t size() { return _dataPoints.size(); }
	const TDataPoint &operator[](const size_t index) const { return _dataPoints[index]; }
	std::vector<TDataPoint>::const_iterator cbegin() const { return _dataPoints.cbegin(); }
	std::vector<TDataPoint>::const_iterator cend() const { return _dataPoints.cend(); }
};

class TSpecies {
private:
	TDataSeries _data;
	std::string _name;
public:
};

#endif /* CORE_TSPECIES_H_ */
