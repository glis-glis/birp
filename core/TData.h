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

class TModelBase;
class TData;
class TDistributionBase;

//-------------------------------------------
// TTimepoints
//-------------------------------------------

class TTimepoints {
private:
	// counts
	std::vector<TypeCounts> _counts_per_species;

	// covariates
	std::vector<TypeCovariateEffort> _covariatesEffort;
	std::vector<TypeCovariateDetection> _covariatesDetection;

	// timepoint id
	size_t _timepoint_id = 0;

	// store effort (updated in MCMC)
	TypeEffort _effort        = 0;
	TypeEffort _effort_old    = 0;
	double _detectionTerm     = 0.0;
	double _detectionTerm_old = 0.0;
	double _effortTerm        = 0.0;
	double _effortTerm_old    = 0.0;

	// covariates
	void _standardizeCovariatesDetection(double MeanDetectionCovariates, double SdDetectionCovariates);
	void _standardizeCovariatesEffort();
	double _calculateEffortTerm(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha);
	double _calculateDetectionTerm(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                               const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);

public:
	TTimepoints() = default;

	TTimepoints(const std::vector<TypeCounts> &Counts, const std::vector<TypeCovariateEffort> &CovsEffort,
	            const std::vector<TypeCovariateDetection> &CovsDetection, const size_t &Timepoint_id);

	explicit TTimepoints(size_t timepoint_id, size_t numSpecies, size_t numCovariatesEffort,
	                     size_t numCovariatesDetection);

	// setters
	void setTimepoint(size_t timepointToAdd);

	void setCovariateEffort(size_t c, TypeCovariateEffort Covariate);
	void setCovariateDetection(size_t c, TypeCovariateDetection Covariate);
	void standardizeCovariates(double MeanDetectionCovariates, double SdDetectionCovariates);

	void write(coretools::TOutputFile &File, std::string location, const TUniqueContainer<size_t> &Timepoints) const;

	// getters
	const std::vector<TypeCovariateEffort> &covariatesEffort() const;
	const std::vector<TypeCovariateDetection> &covariatesDetection() const;
	size_t numSpecies() const;
	TypeCounts counts_per_species(size_t species_id) const;
	TypeEffort effort() const;
	TypeEffort effort_old() const;
	size_t timepoint_id() const;
	bool inferAlpha() const;
	bool inferBeta() const;

	// update effort
	void initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void updateAlphaEffort(TypeAlpha NewAlpha_c1, TypeAlpha OldAlpha_c1, TypeAlpha NewAlpha_c2, TypeAlpha OldAlpha_c2,
	                       size_t c1, size_t c2);
	void updateBetaEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void resetAlphaEffort();
	void resetBetaEffort();

	// functions for simulating covariates
	template<class Object, bool SimulateCovEffort> void simulateCovariate(size_t c, Object obj) {
		if (SimulateCovEffort) {
			_covariatesEffort[c] = obj.sample();
		} else {
			// TODO: how to simulate covariates of detection?
			_covariatesDetection[c] = obj.sample();
		}
	}

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

	TLocations(size_t location_id, size_t TimepointSize, size_t numSpecies, size_t numCovariatesEffort,
	           size_t numCovariatesDetection);

	// setters
	void add(const TTimepoints &pt);

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
	void initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void updateAlphaEffort(TypeAlpha NewAlpha_c1, TypeAlpha OldAlpha_c1, TypeAlpha NewAlpha_c2, TypeAlpha OldAlpha_c2,
	                       size_t c1, size_t c2);
	void updateBetaEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                      const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void resetAlphaEffort();
	void resetBetaEffort();

	// functions for simulating effort
	template<class Object, bool SimulateCovEffort> void simulateCovariate(size_t c, Object obj) {
		for (auto &timepoint : _timepoints) { timepoint.simulateCovariate<Object, SimulateCovEffort>(c, obj); }
	}

	// functions for simulating under Poisson model
	void simulatePoisson(size_t N0, const TModelBase &Model, double detectionProbs);

	// functions for simulating under Negative Binomial model
	void simulateNB(const TModelBase &Model, double a, double mu_ij);
};

//-------------------------------------------
// TMethods
//-------------------------------------------

class TMethods {
private:
	std::string _methodName;
	size_t _method_id;

	std::vector<TLocations> _locations;

	std::vector<size_t> _covariateEffortIDsinUniqueContainer;
	std::vector<size_t> _covariateDetectionIDsinUniqueContainer;

	// variables to identify correct species in epoch model
	std::vector<size_t> _speciesIDsinUniqueContainer;
	std::vector<bool> _hasDataForSpeciesID;
	std::vector<size_t> _indexInCounts;

	// alpha and beta (regression coefficients of covariates)
	std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> _alpha;
	std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> _beta;
	std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> _beta0;

	// update functions for beta
	void _updateBeta(size_t c1, size_t c2, const TData &Data, const std::vector<std::unique_ptr<TModelBase>> &Models);
	void _setValuesForSingleCovariate();

	// simulate
	void _simulateOneCovariateEffort(size_t c, const std::string &Cov);
	void _simulateOneCovariateDetection(size_t c, const std::string &Cov);
	std::vector<std::string> _splitVecWithDistributions(const std::string &Input);
	std::pair<std::string, std::string> _getParamsAndDistr(const std::string &Cov);
	template<class Distribution, bool SimulateCovEffort>
	void _simulateCovariateFromDistribution(size_t c, const std::string &Params) {
		Distribution distr(Params);
		for (size_t i = 0; i < size(); i++) {
			_locations[i].simulateCovariate<Distribution, SimulateCovEffort>(c, distr);
		}
	}
	void _simulateCovariatesDetection();
	void _simulateCovariatesEffort();
	void _simulateCovariates();

public:
	explicit TMethods(const std::string_view &methodName, size_t index,
	                  const std::vector<std::size_t> &speciesIDsinUniqueContainer,
	                  const std::vector<std::size_t> &covariateEffortIDsinUniqueContainer,
	                  const std::vector<std::size_t> &covariateDetectionIDsinUniqueContainer);

	TMethods(const std::string &methodName, size_t locationSize, size_t TimepointSize, size_t index,
	         const std::vector<size_t> &speciesIDsinUniqueContainer,
	         const std::vector<std::size_t> &covariateEffortIDsinUniqueContainer,
	         const std::vector<std::size_t> &covariateDetectionIDsinUniqueContainer);

	TMethods() = default;

	// add parameters alpha and beta
	void addAlphaBeta(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
	                  const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
	                  std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta);
	void initialize(const TUniqueContainer<std::string> &CovariateEffortNames,
	                const TUniqueContainer<std::string> &CovariateDetectionNames);

	// add data
	void add(const TTimepoints &pt, size_t location_id);

	// finalize data filling
	void sorttimes(const std::vector<size_t> &index);
	void standardizeCovariates();
	void fillIDVectors(size_t totalSizeOfUniqueContainer);
	void initializeEffort();

	// getters
	size_t size() const;
	const std::string &name() const;
	size_t method_id() const;
	size_t numCovariatesEffort() const;
	size_t numCovariatesDetection() const;
	bool inferAlpha() const;
	bool inferBeta() const;

	// access counts and covariates
	const std::vector<size_t> &speciesIDsinUniqueContainer() const;
	const std::vector<size_t> &covariateEffortIDsinUniqueContainer() const;
	size_t covariateEffortIDinUniqueContainer(size_t c) const;
	const std::vector<size_t> &covariateDetectionIDsinUniqueContainer() const;
	size_t covariateDetectionIDinUniqueContainer(size_t c) const;
	bool hasDataForSpeciesID(size_t speciesIDinUniqueContainer) const;
	size_t getIndexInCounts(size_t speciesIDinUniqueContainer) const;

	// iterators
	const TLocations &operator[](const size_t index) const { return _locations[index]; }
	std::vector<TLocations>::const_iterator cbegin() const;
	std::vector<TLocations>::iterator begin();
	std::vector<TLocations>::const_iterator cend() const;
	std::vector<TLocations>::iterator end();

	// update effort
	void estimateInitialBeta();
	void updateBeta(const TData &Data, const std::vector<std::unique_ptr<TModelBase>> &Models);
	// functions for simulating
	void simulatePoisson(size_t N0, const TModelBase &Model, const TUniqueContainer<std::string> &locations,
	                     bool SimulateCovariates);

	std::vector<double> simulateNB(size_t N0, const TModelBase &Model, const TUniqueContainer<std::string> &locations,
	                               double A, bool SimulateCovariates);

	// writing
	void write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
	           const TUniqueContainer<size_t> &Timepoints) const;
};

//-------------------------------------------
// TData
//-------------------------------------------

class TData {
private:
	std::vector<TMethods> _methods;

	size_t _numSpecies = 0;

	std::vector<size_t> _numMethLoc; // per species

	// indices (per species)
	std::vector<std::vector<std::pair<size_t, size_t>>> _linear_to_i_j;
	std::vector<std::vector<std::vector<size_t>>> _i_j_to_linear;

	// fill indices
	void _fillNumMethLocs();
	std::pair<size_t, size_t> _getMethAndLocIndex(size_t ij_index, size_t speciesID) const;

public:
	TData();
	virtual ~TData() = default;

	void addMethod(const TMethods &Method);

	void fillMethLocIndices(size_t NumSpecies);

	// iterators
	size_t size() const;
	const TMethods &operator[](const size_t index) const;
	TMethods &operator[](const size_t index);
	std::vector<TMethods>::const_iterator cbegin() const;
	std::vector<TMethods>::iterator begin();
	std::vector<TMethods>::const_iterator cend() const;
	std::vector<TMethods>::iterator end();
	TMethods &back();

	// getters
	size_t numMethLoc(size_t SpeciesIDInUniqueContainer) const;
	size_t i_j_to_linear(size_t SpeciesIDInUniqueContainer, size_t i, size_t j) const;
	const std::pair<size_t, size_t> &linear_to_i_j(size_t SpeciesIDInUniqueContainer, size_t linear) const;
};

#endif // BIRP_TDATA_H
