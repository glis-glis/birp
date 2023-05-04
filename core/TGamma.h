//
// Created by caduffm on 4/17/23.
//

#ifndef BIRP_TGAMMA_H
#define BIRP_TGAMMA_H

#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

class TModelBase;
class TData;

//------------------------------------------
// TGamma
//------------------------------------------

class TGamma {
	// wraps parameters gamma and timesOfChange
protected:
	// stattools parameters
	stattools::TParameterTyped<TypeGamma, 1> *_gamma        = nullptr;
	stattools::TParameterTyped<TypeTime, 1> *_timesOfChange = nullptr;

	// number of epochs
	size_t _numEpochs = 0;
	std::vector<size_t> _epochStarts;

	// rho
	coretools::TMatrix<size_t> _rho;

	// Jeffrey's prior
	double _tryJP = 0.0;
	double _curJP = 0.0;

	// initialization
	void _initializeEpochStarts(const TUniqueContainer<size_t> &Timepoints);
	void _fillRho(const TUniqueContainer<size_t> &Timepoints);
	bool _timePointIsWithinEpoch(const TUniqueContainer<size_t> &Timepoints, size_t t, size_t m) const;
	auto _get_Y_m(const TData &Data, const TUniqueContainer<size_t> &Timepoints, size_t m, size_t SpeciesID) const;
	arma::vec _doOLS(const std::vector<double> &Y_m, const std::vector<TypeTime> &T_m,
	                 const std::vector<size_t> &DataSeriesIndex, size_t L) const;
	void _augmentYWithOLSEstimatesStart(size_t curEpoch, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                    std::vector<std::vector<TypeTime>> &t,
	                                    std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const;
	void _augmentYWithOLSEstimatesEnd(size_t curEpoch, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                  std::vector<std::vector<TypeTime>> &t,
	                                  std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const;
	void _augmentDataSeriesIndexWithOLSEstimates(std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij, size_t l,
	                                             std::vector<std::vector<size_t>> &dataSeriesIndex, size_t curEpoch,
	                                             size_t otherEpoch) const;
	bool _dontAugmentDataSeries(const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                            const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                            const std::vector<std::vector<TypeTime>> &t, size_t l, size_t curEpoch,
	                            size_t otherEpoch, TypeTime EpochBorder) const;
	bool _epochIsIdentifiable(size_t L, const std::vector<TypeTime> &t,
	                          const std::vector<size_t> &dataSeriesIndex) const;

public:
	TGamma(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TGamma() = default;

	void initialize(const TUniqueContainer<size_t> &Timepoints);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<size_t> &Timepoints, size_t SpeciesID);
	void initializeJP(const TData &Data, const TModelBase &Model);

	std::vector<size_t> getEpochStarts() const;

	double getLogPriorRatio(const TData &Data, const TModelBase &Model);

	bool updateSpecificIndex(size_t m);
	bool acceptOrReject(double logH, size_t m);

	void setSimulatedGammaToZero();

	size_t numEpochs() const;
	size_t numTimepoints() const;
	TypeGamma gamma(size_t m) const;
	size_t rho(size_t row, size_t col) const;
	std::vector<double> phi() const;
};

//------------------------------------------
// TStochasticPrior
//------------------------------------------

class TStochasticPrior {
	// yellow box
protected:
	stattools::TParameterTyped<TypeLogPhi, 2> *_logPhi = nullptr;
	stattools::TParameterTyped<TypeSigma, 1> *_sigma   = nullptr;

	// TGamma
	TGamma _gamma;

	// rho
	coretools::TMatrix<size_t> _stochasticRho;
	std::vector<size_t> _deltaTimeVec;

	// temporary variable: sum of rho * gamma
	std::vector<double> _newSumRhoGamma;
	std::vector<double> _oldSumRhoGamma;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<double>> _tryLogTransitionProbs;
	std::vector<std::vector<double>> _curLogTransitionProbs;

	// fill rho and delta
	void _fillStochasticRho(const TUniqueContainer<size_t> &Timepoints);
	void _initializeSumRhoGamma();
	void _initializeTransitionProbabilities();

	// transition probabilities
	double _calculateLogTransitionProbability(size_t j, size_t k) const;
	coretools::probdist::TNormalDistr _getMeanSdTransitionProbability(size_t j, size_t k) const;

public:
	TStochasticPrior(stattools::TParameterTyped<TypeLogPhi, 2> *LogPhi, stattools::TParameterTyped<TypeSigma, 1> *Sigma,
	                 stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	                 stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TStochasticPrior() = default;

	void initialize(const TUniqueContainer<size_t> &Timepoints);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<size_t> &Timepoints, size_t SpeciesID);
	void initializeTempVariables(const TData &Data, const TModelBase &Model, bool CalculateJP);

	// update gamma
	void updateGamma(const TData &Data, const TModelBase &Model);
	double calculateLL() const;
	double calculateLLRatioForUpdateGamma();
	void swapTryCur();
	void swapTryCur(size_t j, size_t k);

	double calculateLogPriorRatio(size_t j, size_t k);

	// simulations
	void setSimulatedGammaToZero();
	double sampleFromTransitionProbability(size_t j, size_t k) const;

	// getters
	std::vector<double> deterministicPhi() const;
};

#endif // BIRP_TGAMMA_H
