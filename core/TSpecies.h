//
// Created by liam on 12/09/22.
//

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include "BirpTypes.h"
#include <algorithm>
#include <vector>
#include "coretools/algorithms.h"

template<typename T> class TUniqueContainer {
protected:
	std::vector<T> _names;

public:
	bool exists(const T &Value) const { return std::find(_names.begin(), _names.end(), Value) != _names.end(); }

	size_t add(const T &Value) {
		if (!exists(Value)) {
			_names.push_back(Value);
			return _names.size() - 1;
		}
		return getIndex(Value);
	}

	size_t getIndex(const T &Value) const {
		auto it = std::find(_names.begin(), _names.end(), Value);
		if (it != _names.end()) { return it - _names.begin(); }
		DEVERROR("Value ", Value, " not found in vector. This should not happen.");
	}

	size_t size() const { return _names.size(); }
	void resize(size_t size) { _names.resize(size); }

	const T &operator[](size_t Index) const { return _names[Index]; }

	typename std::vector<T>::const_iterator cbegin() const { return _names.cbegin(); }
	typename std::vector<T>::iterator begin() { return _names.begin(); }
	typename std::vector<T>::const_iterator cend() const { return _names.cend(); }
	typename std::vector<T>::iterator end() { return _names.end(); }

	bool empty() const { return _names.size() == 0; }
	T back() const { return _names.back(); }
};

class TDataPoint {
public:
	size_t timepoint_id = 0;
	TypeCounts counts   = 0;
	TypeEffort effort   = 0;

	TDataPoint(size_t Timpoint_id, TypeCounts Counts, TypeEffort Effort) {
		timepoint_id = Timpoint_id;
		counts       = Counts;
		effort       = Effort;
	}
	TDataPoint() = default;

	bool operator<(const TDataPoint &other) const { return timepoint_id < other.timepoint_id; }
};

class TDataSeries {
private:
	std::vector<TDataPoint> _dataPoints;
	size_t _location_id = 0;
	size_t _method_id   = 0;
	TypeCounts _nu; // sums of counts at one location and method (ij) across all time points (k)

public:
	TDataSeries(size_t Location_id, size_t Method_id) {
		_location_id = Location_id;
		_method_id   = Method_id;
	}

	TDataSeries(const TDataPoint &Datapoint, size_t Location_id, size_t Method_id) {
		_location_id = Location_id;
		_method_id   = Method_id;
		add(Datapoint);
	}

	// set
	void add(const TDataPoint &Datapoint) {
		_dataPoints.push_back(Datapoint);
		_nu += Datapoint.counts;
	}

	// get
	size_t size() const { return _dataPoints.size(); }

	const TDataPoint &operator[](const size_t index) const { return _dataPoints[index]; }

	std::vector<TDataPoint>::const_iterator cbegin() const { return _dataPoints.cbegin(); }
	std::vector<TDataPoint>::iterator begin() { return _dataPoints.begin(); }
	std::vector<TDataPoint>::const_iterator cend() const { return _dataPoints.cend(); }
	std::vector<TDataPoint>::iterator end() { return _dataPoints.end(); }

	// getters
	size_t location_id() const { return _location_id; }
	size_t method_id() const { return _method_id; }
	TypeCounts nu() const {return _nu;}

	void sorttimes(const std::vector<size_t> &index) {
		std::vector<size_t> positionOfTimepointIds;
		for (size_t i = 0; i < _dataPoints.size(); i++) {
			auto it    = std::find(index.begin(), index.end(), _dataPoints[i].timepoint_id);
			size_t pos = it - index.begin();
			positionOfTimepointIds.push_back(pos);
		}
		// sort the newly found positions of timepoints using ranksort()
		std::vector<size_t> sortingIndex = coretools::rankSort(positionOfTimepointIds);
		// sort _dataPoints
		_dataPoints                      = coretools::sortContainerByRank(_dataPoints, sortingIndex);
		// now sort positionOfTimepointIds
		positionOfTimepointIds           = coretools::sortContainerByRank(positionOfTimepointIds, sortingIndex);
		// now set timepoint_id
		for (size_t i = 0; i < _dataPoints.size(); i++) { _dataPoints[i].timepoint_id = positionOfTimepointIds[i]; }
        // issue warning if any DataSeries contains less than two time points with nonzero efforts
        size_t nonZeroEffort = 0;
        for (size_t i = 0; i < _dataPoints.size(); i++) {
            if(_dataPoints[i].effort > 0){ nonZeroEffort += 1; }
        }
        if(_dataPoints.size() < 2 || nonZeroEffort < 2){
            // TODO write to log file
            // TLogFile.write("What a noob")
        }
	}
};

class TSpecies {
private:
	std::vector<TDataSeries> _data;
	std::string _name;

public:
	explicit TSpecies(const std::string &Name) { _name = Name; }

	explicit TSpecies(const std::string &Name, const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex) {
		_name = Name;
		add(pt, LocationIndex, MethodIndex);
	}

	std::string name() const { return _name; }

	void add(const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex) {
		bool dataPointExists = false;
		for (size_t d = 0; d < _data.size(); d++) {
			if (_data[d].location_id() == LocationIndex && _data[d].method_id() == MethodIndex) {
				dataPointExists = true;
				_data[d].add(pt);
				break;
			}
		}
		if (!dataPointExists) { _data.emplace_back(pt, LocationIndex, MethodIndex); }
	}


	void sorttimes(const std::vector<size_t> &index) {
		for (size_t ds = 0; ds < _data.size(); ds++) { _data[ds].sorttimes(index); }
	}

	size_t size() const {return _data.size();}

	std::vector<TDataSeries>::const_iterator cbegin() const { return _data.cbegin(); }
	std::vector<TDataSeries>::iterator begin() { return _data.begin(); }
	std::vector<TDataSeries>::const_iterator cend() const { return _data.cend(); }
	std::vector<TDataSeries>::iterator end() { return _data.end(); }

	const TDataSeries &operator[](const size_t index) const { return _data[index]; }
};

#endif /* CORE_TSPECIES_H_ */
