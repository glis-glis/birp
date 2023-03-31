//
// Created by liam on 10.03.23.
//

#ifndef BIRP_TDATA_H
#define BIRP_TDATA_H

#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <string>
#include <vector>

//-------------------------------------------
// TTimepoints
//-------------------------------------------

class TTimepoints {
private:
	std::vector<TypeCovariate> _covariates;
	std::vector<TypeCounts> _counts_per_species;
	TypeEffort _effort     = 0; // TODO how to store effort and distinguish betewen covariate?
	TypeEffort _effort_old = 0;
	size_t _timepoint_id   = 0;

	void _simulatePoisson(double lambda);

public:
	TTimepoints() = default;

	TTimepoints(const std::vector<TypeCounts> &countVec, const std::vector<TypeCovariate> &covVec,
	            const size_t &timepoint_id);

	explicit TTimepoints(size_t timepoint_id, size_t numSpecies, size_t numCovariates);

	// setters
	void setTimepoint(size_t timepointToAdd);

	void setEffort(TypeEffort effortToAdd);

	void write(coretools::TOutputFile &File, std::string location, const TUniqueContainer<size_t> &Timepoints) const;

	// getters
	std::vector<TypeCovariate> covariates();

	size_t numSpecies() const;

	TypeCounts counts_per_species(size_t species_id) const;

	std::vector<size_t> species_id();

	TypeEffort effort() const;

	TypeEffort effort_old();

	size_t timepoint_id() const;

	// functions for simulating effort
	template<class Object> void simulateEffort(Object obj) {
		_effort = obj.sample();
		if (_covariates.size() == 1) { _covariates[0] = _effort; } // TODO is this correct
	}

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, double phi, double detectionProbs);

	void simulatePoisson(double phi);

	// functions for simulating under Negative Binomial model
	void simulateNB(double phi, double a, double mu_ij);
};

//-------------------------------------------
// TLocations
//-------------------------------------------

class TLocations {
private:
	size_t _location_id = 0;
	std::vector<TTimepoints> _timepoints;
	std::vector<TypeNu> _nu;

public:
	TLocations() = default;

	TLocations(const TTimepoints &pt, size_t location_id);

	TLocations(size_t location_id, size_t TimepointSize, size_t numSpecies, size_t numCovariates);

	// setters
	void add(const TTimepoints &pt);

	void checkEffortBiggerThanZero();

	void sorttimes(const std::vector<size_t> &index);

	// getters
	size_t size() const;

	TypeNu nu(size_t species_id) const;

	size_t location_id() const;

	const TTimepoints &operator[](const size_t index) const;

	void write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<size_t> &Timepoints) const;

	std::vector<TTimepoints>::const_iterator cbegin() const;

	std::vector<TTimepoints>::iterator begin();

	std::vector<TTimepoints>::const_iterator cend() const;

	std::vector<TTimepoints>::iterator end();

	// functions for simulating effort
	template<class Object> void simulateEffort(Object obj) {
		for (auto &timepoint : _timepoints) { timepoint.simulateEffort(obj); }
	}

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, const std::vector<double> &phi, double detectionProbs);

	void simulatePoisson(const std::vector<double> &phi);

	// functions for simulating under Negative Binomial model
	void simulateNB(const std::vector<double> &phi, double a, double mu_ij);
};

//-------------------------------------------
// TMethods
//-------------------------------------------

class TMethods {
private:
	std::string _methodName;
	std::vector<double> _beta;
	std::vector<TLocations> _locations;
	size_t _method_id;
	std::vector<std::string> _covariateNames;
	// variables to identify correct species in epoch model
	std::vector<size_t> _speciesIDsinUniqueContainer;
	std::vector<bool> _hasDataForSpeciesID;
	std::vector<size_t> _indexInCounts;

	template<class Distribution> void _simulateEffort(const std::string &Params) {
		Distribution distr(Params);
		for (size_t i = 0; i < size(); i++) { _locations[i].simulateEffort(distr); }
	}

	void _simulateEffort();
	std::pair<size_t, size_t>
	_getIndexInDetectionProbabilitiesToSimulate(const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	                                            const TUniqueContainer<std::string> &locations, size_t i);

public:
	explicit TMethods(const std::string_view &methodName, size_t index, const std::vector<std::size_t> &species_id);

	TMethods(const std::string &methodName, size_t locationSize, size_t TimepointSize, size_t index,
	         const std::vector<size_t> &species_ids, const std::vector<std::string> &covariateNames);

	TMethods() = default;

	void fillIDVectors(size_t totalSizeOfUniqueContainer);

	// getters
	size_t size() const;
	std::string methodName() const;
	std::vector<std::string> covariateNames() const;
	size_t method_id() const;
	const std::vector<size_t> &speciesIDsinUniqueContainer() const;
	bool hasDataForSpeciesID(size_t speciesIDinUniqueContainer) const;
	size_t getIndexInCounts(size_t speciesIDinUniqueContainer) const;

	void add(const TTimepoints &pt, size_t location_id);

	std::string name() const;

	void sorttimes(const std::vector<size_t> &index);

	void checkEffortBiggerThanZero();

	void write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<size_t> &Timepoints) const;

	const TLocations &operator[](const size_t index) const { return _locations[index]; }

	std::vector<TLocations>::const_iterator cbegin() const;

	std::vector<TLocations>::iterator begin();

	std::vector<TLocations>::const_iterator cend() const;

	std::vector<TLocations>::iterator end();

	// functions for simulating
	void simulatePoisson(size_t N0, const std::vector<double> &phi,
	                     const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	                     const TUniqueContainer<std::string> &locations);

	void simulatePoisson(const std::vector<double> &phi);

	std::vector<double> simulateNB(size_t N0, const std::vector<double> &phi,
	                               const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	                               const TUniqueContainer<std::string> &locations, double A);

	void simulateNB(const std::vector<double> &phi, double A, const std::vector<double> &Mu, bool SimulateEffort);
};

#endif // BIRP_TDATA_H
