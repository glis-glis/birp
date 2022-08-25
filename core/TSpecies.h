/*
 * TSpecies.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include <algorithm>
#include <vector>
#include "BirpTypes.h"


struct TDataPoint{
	TypeIndex location;
	TypeIndex method;
	TypeIndex timepoint;
	TypeCounts counts;
	TypeEffort effort;

	bool operator<(const TDataPoint & other){
		if(location < other.location){
			return true;
		} else if(location > other.location){
			return false;
		} else if(method < other.method){
			return true;
		} else if(method > other.method){
			return false;
		} else if(timepoint < other.timepoint){
			return true;
		} else {
			return false;
		}
	}
};

class TData{
private:
	std::vector<TDataPoint> _dataPoints;

public:
	TData() = default;

	TData(const std::string Filename){
		//read
		//TODO

		//ensure data is sorted
		std::sort(_dataPoints.begin(), _dataPoints.end());

		//check if some location and methods have only one time point -> remove them!
	}

	size_t size(){ return _dataPoints.size(); }

	const TDataPoint& operator[](const size_t index) const {
		return _dataPoints[index];
	}

	std::const_iterator< std::vector<TDataPoint> > cbegin() const {
		return _dataPoints.cbegin();
	}

	std::const_iterator< std::vector<TDataPoint> > cend() const {
		return _dataPoints.cend();
	}
};

class TSpecies{
private:
	TData _data;
public:


};


#endif /* CORE_TSPECIES_H_ */
