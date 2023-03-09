//
// Created by liam on 12/09/22.
//

#ifndef CORE_TSPECIES_H_
#define CORE_TSPECIES_H_

#include "BirpTypes.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include "coretools/Storage/TStorage.h"
#include <algorithm>
#include <vector>
#include "TUniqueContainer.h"
#include "TCovariates.h"

//-------------------------------------------
// TZeroPoisDistr
//-------------------------------------------

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

	size_t sample(coretools::Probability p, coretools::StrictlyPositive<double> lambda) const;

	// functions using member parameters
	size_t sample() const noexcept;
};

//-------------------------------------------
// TDataPoint
//-------------------------------------------

class TDataPoint {
private:
	void _simulatePoisson(double lambda);

public:
	size_t timepoint_id = 0;
	TypeCounts counts   = 0;
	TypeEffort effort   = 0;
    TypeEffort effort_old   = 0;

	TDataPoint(size_t Timpoint_id, TypeCounts Counts, TypeEffort Effort);
	TDataPoint() = default;
	TDataPoint(size_t Timpoint_id);
	TDataPoint (size_t Timepoint_id, TypeCounts Counts);

	// functions for simulating effort
	template<class Object> void simulateEffort(Object obj) { effort = obj.sample(); }

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, double phi, double detectionProbs);
	void simulatePoisson(double phi);

	// functions for simulating under Negative Binomial model
	void simulateNB(double phi, double a, double mu_ij);

	// write to file
	void write(coretools::TOutputFile &File, const std::string &SpeciesName, const std::string &Location,
	           const std::string &Method, const TUniqueContainer<size_t> &Timepoints) const;

	bool operator<(const TDataPoint &other) const;

    bool timepointExists(size_t timepointIndex);
};

//-------------------------------------------
// TDataSeries
//-------------------------------------------

class TDataSeries {
private:
	std::vector<TDataPoint> _dataPoints;
	size_t _location_id = 0;
	size_t _method_id   = 0;
	TypeNu _nu; // sums of counts at one location and method (ij) across all time points (k)

public:
	explicit TDataSeries(size_t Location_id, size_t Method_id);
	explicit TDataSeries(const TDataPoint &Datapoint, size_t Location_id, size_t Method_id);
	explicit TDataSeries(size_t Location_id, size_t Method_id, size_t numTimepoints);

	void sorttimes(const std::vector<size_t> &index);

	void checkData(const std::string &name, const TUniqueContainer<std::string> &locations,
	               const TUniqueContainer<std::string> &methods);

    bool IJK_exists(size_t methodIndex, size_t locationIndex, size_t timepointIndex);
    void checkCovExists(TCovariate covariate);

	// setters
	void add(const TDataPoint &Datapoint);

	// getters
	size_t size() const;

	size_t location_id() const;
	size_t method_id() const;
	TypeNu nu() const;

	const TDataPoint &operator[](const size_t index) const { return _dataPoints[index]; }

	std::vector<TDataPoint>::const_iterator cbegin() const;
	std::vector<TDataPoint>::iterator begin();
	std::vector<TDataPoint>::const_iterator cend() const;
	std::vector<TDataPoint>::iterator end();

	// functions for simulating effort
	template<class Object> void simulateEffort(Object obj) {
		for (auto &dataPoint : _dataPoints) { dataPoint.simulateEffort(obj); }
	}

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, const std::vector<double> &phi, double detectionProbs);
	void simulatePoisson(const std::vector<double> &phi);

	// functions for simulating under Negative Binomial model
	void simulateNB(const std::vector<double> &phi, double a, double mu_ij);

	// write to file
	void write(coretools::TOutputFile &File, const std::string &SpeciesName,
	           const TUniqueContainer<std::string> &Locations, const TUniqueContainer<std::string> &Methods,
	           const TUniqueContainer<size_t> &Timepoints) const;
};

//-------------------------------------------
// TSpecies
//-------------------------------------------

class TSpecies {
private:
	std::vector<TDataSeries> _data;
	std::string _name;

	template<class Distribution> void _simulateEffort(const std::string &Params) {
		Distribution distr(Params);
		for (size_t i = 0; i < size(); i++) { _data[i].simulateEffort(distr); }
	}

	void _simulateEffort();

public:
	explicit TSpecies(const std::string &Name);
	explicit TSpecies(const std::string &Name, const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex);
	explicit TSpecies(const std::string &Name, size_t numLocations, size_t numMethods, size_t numTimepoints);

	// setters
	void add(const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex);

	void checkData(const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods);
	void sorttimes(const std::vector<size_t> &index);
    bool IJK_exists(size_t methodIndex, size_t locationIndex, size_t timepointIndex);
    void checkCovExists(const TCovariate &covariate);
    void setFixedEffort(const TCovariate &covariate);

	// getters
	std::string name() const;
	size_t size() const;

	const TDataSeries &operator[](const size_t index) const { return _data[index]; }

	std::vector<TDataSeries>::const_iterator cbegin() const;
	std::vector<TDataSeries>::iterator begin();
	std::vector<TDataSeries>::const_iterator cend() const;
	std::vector<TDataSeries>::iterator end();

	// functions for simulating
	std::pair<size_t,size_t> getIndexToSimulate(const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	                                  const TUniqueContainer<std::string> &locations,
	                                  const TUniqueContainer<std::string> &methods, size_t i);
	void simulatePoisson(size_t N0, const std::vector<double> &phi,
	               const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	               const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods);

	void simulatePoisson(const std::vector<double> &phi);
	void simulateNB(const std::vector<double> &phi, double A, const std::vector<double> &Mu, bool SimulateEffort);
	std::vector<double> simulateNB(size_t N0, const std::vector<double> &phi, const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	                          const TUniqueContainer<std::string> &locations,
	                          const TUniqueContainer<std::string> &methods, double A);

	void write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<std::string> &Methods, const TUniqueContainer<size_t> &Timepoints) const;
};

#endif /* CORE_TSPECIES_H_ */
