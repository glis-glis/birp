//
// Created by liam on 12/09/22.
//

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include "BirpTypes.h"
#include <algorithm>
#include <vector>

template<typename T>
class TUniqueContainer{
protected:
	std::vector<T> _names;
public:
	bool exists(const T &Value) const {
		return std::find(_names.begin(), _names.end(),Value)!=_names.end();
	}
	void add(const T &Value){
		if (!exists(Value)) {
			_names.push_back(Value);
		}
	}
	size_t getIndex(const T &Value) const {
		auto it = std::find(_names.begin(), _names.end(), Value);
		if (it != _names.end()) {
			return  it - _names.begin();
		}
		DEVERROR("Value ",Value, " not found in vector. This should not happen.");
	}

	size_t size() const {
		return _names.size();
	}

	const T &operator[](size_t Index) const { return _names[Index]; }
	typename std::vector<T>::const_iterator cbegin() const { return _names.cbegin(); }
	typename std::vector<T>::const_iterator cend() const { return _names.cend(); }

};

class TDataPoint {
public:
	size_t timepoint_id = 0;
	TypeCounts counts = 0;
	TypeEffort effort = 0;
    TDataPoint(size_t Timpoint_id, TypeCounts Counts, TypeEffort Effort){
		timepoint_id = Timpoint_id;
		counts = Counts;
		effort = Effort;
	}
	bool operator<(const TDataPoint &other) const { return timepoint_id < other.timepoint_id; }
};

class TDataSeries {
private:
	std::vector<TDataPoint> _dataPoints;
	size_t _location_id = 0;
	size_t _method_id = 0;

public:
	TDataSeries(size_t Location_id, size_t Method_id){
		_location_id = Location_id;
		_method_id = Method_id;
	}

    // set
    void add(const TDataPoint &Datapoint){
		_dataPoints.push_back(Datapoint);
	}


	// get
	size_t size() { return _dataPoints.size(); }
	const TDataPoint &operator[](const size_t index) const { return _dataPoints[index]; }
	std::vector<TDataPoint>::const_iterator cbegin() const { return _dataPoints.cbegin(); }
	std::vector<TDataPoint>::const_iterator cend() const { return _dataPoints.cend(); }


};

class TSpecies {
private:
	std::vector<TDataSeries> _data;
public:

};



#endif /* CORE_TSPECIES_H_ */
