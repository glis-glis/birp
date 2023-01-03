//
// Created by liam on 12/09/22.
//

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include "BirpTypes.h"
#include "coretools/algorithms.h"
#include <algorithm>
#include <vector>
#include "coretools/Strings/splitters.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include "coretools/Storage/TStorage.h"

template<typename T> class TUniqueContainer {
protected:
	std::vector<T> _names;

public:
	bool exists(const T &Value) const {
		return std::find(_names.begin(), _names.end(), Value) != _names.end();
	}

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

// ------------------------------------------------------ //

class TDataPoint {
private:
	void _simulatePoisson(double lambda);

public:
	size_t timepoint_id = 0;
	TypeCounts counts   = 0;
	TypeEffort effort   = 0;

	TDataPoint(size_t Timpoint_id, TypeCounts Counts, TypeEffort Effort);
	TDataPoint() = default;
	TDataPoint(size_t Timpoint_id);

	void simulate(size_t N0, double phi, double detectionProbs);
	template<class Object> void simulate(Object obj, double phi){
		effort = obj.sample();
		_simulatePoisson(phi);
	}
	void simulate(double effort, double phi);

	bool operator<(const TDataPoint &other) const;
};

// ------------------------------------------------------ //

class TDataSeries {
private:
	std::vector<TDataPoint> _dataPoints;
	size_t _location_id = 0;
	size_t _method_id   = 0;
	TypeCounts _nu; // sums of counts at one location and method (ij) across all time points (k)

public:
	TDataSeries(size_t Location_id, size_t Method_id);
	TDataSeries(const TDataPoint &Datapoint, size_t Location_id, size_t Method_id);
	TDataSeries(size_t Location_id, size_t Method_id, size_t numTimepoints);
	// set
	void add(const TDataPoint &Datapoint);

	// get
	size_t size() const;

	const TDataPoint &operator[](const size_t index) const { return _dataPoints[index]; }

	std::vector<TDataPoint>::const_iterator cbegin() const;
	std::vector<TDataPoint>::iterator begin();
	std::vector<TDataPoint>::const_iterator cend() const;
	std::vector<TDataPoint>::iterator end();

	// getters
	size_t location_id() const;
	size_t method_id() const;
	TypeCounts nu() const;

	void sorttimes(const std::vector<size_t> &index);

	void checkData(const std::string &name, const TUniqueContainer<std::string> &locations,
	               const TUniqueContainer<std::string> &methods);

	void simulate(size_t N0, const std::vector<double>& phi, double detectionProbs);
	template<class Object> void simulate(Object obj, const std::vector<double>& phi){
		for (size_t t = 0; t < _dataPoints.size(); t++) { _dataPoints[t].simulate(obj, phi[t]); }
	}
	void simulate(double effort, const std::vector<double>& phi);
};

// ------------------------------------------------------ //

class TSpecies {
private:
	std::vector<TDataSeries> _data;
	std::string _name;

public:
	explicit TSpecies(const std::string &Name);
	explicit TSpecies(const std::string &Name, const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex);
	TSpecies(const std::string &Name, size_t numLocations, size_t numMethods, size_t numTimepoints);

	std::string name() const;

	void add(const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex);

	void sorttimes(const std::vector<size_t> &index);

	size_t size() const;

	void checkData(const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods);

	std::vector<TDataSeries>::const_iterator cbegin() const;
	std::vector<TDataSeries>::iterator begin();
	std::vector<TDataSeries>::const_iterator cend() const;
	std::vector<TDataSeries>::iterator end();

	const TDataSeries &operator[](const size_t index) const { return _data[index]; }

	void simulate(size_t N0, const std::vector<double> &phi,
	              const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	              const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods);

	void simulate(const std::vector<double>& phi);


};


//------------------------------------------------
// zero-Poisson distribution
//------------------------------------------------
class TZeroPoisDistr {
private:
	coretools::Probability _p;
	coretools::StrictlyPositive<double> _lambda;

public:
	TZeroPoisDistr(coretools::Probability p, coretools::StrictlyPositive<double> lambda);
	TZeroPoisDistr(std::string_view parameterString) { set(parameterString); }

	void set(coretools::Probability p, coretools::StrictlyPositive<double> Lambda);
	void set(std::string_view parameterString);

	static constexpr std::string_view name = "zeroPoisson";

	static size_t sample( coretools::Probability p, coretools::StrictlyPositive<double> lambda);

	// functions using member parameters
	[[nodiscard]] size_t sample() const noexcept;


};


#endif /* CORE_TSPECIES_H_ */
