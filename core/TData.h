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

class TEpochModel;

//-------------------------------------------
// TTimepoints
//-------------------------------------------

class TTimepoints {
private:
	std::vector<TypeCovariate> _covariates;
	std::vector<TypeCounts> _counts_per_species;
	TypeEffort _effort     = 0;
	TypeEffort _effort_old = 0;
	size_t _timepoint_id   = 0;

public:
	TTimepoints() = default;

	TTimepoints(const std::vector<TypeCounts> &countVec, const std::vector<TypeCovariate> &covVec,
	            const size_t &timepoint_id);

	explicit TTimepoints(size_t timepoint_id, size_t numSpecies, size_t numCovariates);

	// setters
	void setTimepoint(size_t timepointToAdd);

	void setCovariate(size_t c, TypeCovariate Covariate);
	void standardizeCovariatesAndAddIntercept(TypeCovariate Min, TypeCovariate Max);

	void write(coretools::TOutputFile &File, std::string location, const TUniqueContainer<size_t> &Timepoints) const;

	// getters
	std::vector<TypeCovariate> covariates() const;
	size_t numSpecies() const;
	TypeCounts counts_per_species(size_t species_id) const;
	TypeEffort effort() const;
	TypeEffort effort_old() const;
	size_t timepoint_id() const;
	bool inferEffort() const;

	// update effort
	void initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void updateEffort(TypeBeta NewBeta_c1, TypeBeta OldBeta_c1, TypeBeta NewBeta_c2, TypeBeta OldBeta_c2, size_t c1,
	                  size_t c2);
	void resetEffort();

	// functions for simulating covariates
	template<class Object> void simulateCovariate(size_t c, Object obj) { _covariates[c] = obj.sample(); }

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, double phi, double detectionProbs);

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

	// update effort
	void initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void updateEffort(TypeBeta NewBeta_c1, TypeBeta OldBeta_c1, TypeBeta NewBeta_c2, TypeBeta OldBeta_c2, size_t c1,
	                  size_t c2);
	void resetEffort();

	// functions for simulating effort
	template<class Object> void simulateCovariate(size_t c, Object obj) {
		for (auto &timepoint : _timepoints) { timepoint.simulateCovariate(c, obj); }
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
	size_t _method_id;

	std::vector<TLocations> _locations;

	std::vector<size_t> _covariateIDsinUniqueContainer;
	// variables to identify correct species in epoch model
	std::vector<size_t> _speciesIDsinUniqueContainer;
	std::vector<bool> _hasDataForSpeciesID;
	std::vector<size_t> _indexInCounts;

	// beta (regression coefficients)
	std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> _beta;

	// update functions for beta
	void _updateBeta(size_t c1, size_t c2, const std::vector<std::unique_ptr<TEpochModel>> &EpochModels);

	// standardize covariates
	std::pair<TypeCovariate, TypeCovariate> _getMinMaxCovariates() const;
	void _setBetaSingleCovariate();

	// simulate
	void _simulateOneCovariate(size_t c, const std::string &Cov);
	std::vector<std::string> _splitVecWithDistributions(const std::string &Input);
	template<class Distribution> void _simulateCovariateFromDistribution(size_t c, const std::string &Params) {
		Distribution distr(Params);
		for (size_t i = 0; i < size(); i++) { _locations[i].simulateCovariate(c, distr); }
	}
	void _simulateCovariate();
	coretools::Probability _getDetectionProbability(std::tuple<std::vector<std::string>, std::vector<std::string>,
	                                                           std::vector<coretools::Probability>> &detectionProbs,
	                                                const TUniqueContainer<std::string> &locations, size_t j);

public:
	explicit TMethods(const std::string_view &methodName, size_t index,
	                  const std::vector<std::size_t> &speciesIDsinUniqueContainer,
	                  const std::vector<std::size_t> &covariateIDsinUniqueContainer);

	TMethods(const std::string &methodName, size_t locationSize, size_t TimepointSize, size_t index,
	         const std::vector<size_t> &speciesIDsinUniqueContainer,
	         const std::vector<std::size_t> &covariateIDsinUniqueContainer);

	TMethods() = default;

	// add parameters beta
	void addBeta(std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);

	// add data
	void add(const TTimepoints &pt, size_t location_id);

	// finalize data filling
	void sorttimes(const std::vector<size_t> &index);
	void checkEffortBiggerThanZero();
	void standardizeCovariatesAndAddIntercept();
	void fillIDVectors(size_t totalSizeOfUniqueContainer);
	void initializeEffort();

	// getters
	size_t size() const;
	const std::string &name() const;
	size_t method_id() const;
	size_t numCovariates() const;
	bool inferEffort() const;

	// access counts and covariates
	const std::vector<size_t> &speciesIDsinUniqueContainer() const;
	const std::vector<size_t> &covariateIDsinUniqueContainer() const;
	size_t covariateIDinUniqueContainer(size_t c) const;
	bool hasDataForSpeciesID(size_t speciesIDinUniqueContainer) const;
	size_t getIndexInCounts(size_t speciesIDinUniqueContainer) const;

	// iterators
	const TLocations &operator[](const size_t index) const { return _locations[index]; }
	std::vector<TLocations>::const_iterator cbegin() const;
	std::vector<TLocations>::iterator begin();
	std::vector<TLocations>::const_iterator cend() const;
	std::vector<TLocations>::iterator end();

	// update effort
	void initializeInferred(const TUniqueContainer<std::string> &CovariateNames);
	void estimateInitialBeta();
	void updateBeta(const std::vector<std::unique_ptr<TEpochModel>> &EpochModels);

	// functions for simulating
	void simulatePoisson(size_t N0, const std::vector<double> &phi,
	                     std::tuple<std::vector<std::string>, std::vector<std::string>,
	                                std::vector<coretools::Probability>> &detectionProbs,
	                     const TUniqueContainer<std::string> &locations, bool SimulateCovariates);

	std::vector<double> simulateNB(size_t N0, const std::vector<double> &phi,
	                               std::tuple<std::vector<std::string>, std::vector<std::string>,
	                                          std::vector<coretools::Probability>> &detectionProbs,
	                               const TUniqueContainer<std::string> &locations, double A, bool SimulateCovariates);

	// writing
	void write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<size_t> &Timepoints) const;
};

#endif // BIRP_TDATA_H
