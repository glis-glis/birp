#ifndef BIRP_TGAMMA_H
#define BIRP_TGAMMA_H

#include "BirpDAG.h"
#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "coretools/Math/TMatrix.h"
#include <armadillo>
#include <vector>

class TModelBase;
class TData;

//------------------------------------------
// TOLSGamma
//------------------------------------------

class TOLSGamma {
	// class used to initialize gamma with OLS estimate
private:
	// values of gamma
	std::vector<double> _gamma;

	// number of epochs
	size_t _numEpochs = 0;
	std::vector<TypeTime> _epochStarts;

	bool _timePointIsWithinEpoch(const TUniqueContainer<TypeTime> &Timepoints, size_t t, size_t m) const;
	size_t _getNextEpochInitialGamma(const std::vector<bool> &EpochInitialized,
	                                 const std::vector<std::vector<double>> &y,
	                                 const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const;
	auto _get_Y_m(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t m, size_t SpeciesID) const;
	static arma::vec _doOLS(const std::vector<double> &Y_m, const std::vector<TypeTime> &T_m,
	                        const std::vector<size_t> &DataSeriesIndex, size_t L);
	void _augmentYWithOLSEstimatesStart(size_t curEpoch, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                    std::vector<std::vector<TypeTime>> &t,
	                                    std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const;
	void _augmentYWithOLSEstimatesEnd(size_t curEpoch, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                  std::vector<std::vector<TypeTime>> &t,
	                                  std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const;
	static void _augmentDataSeriesIndexWithOLSEstimates(std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                                    size_t l, std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                                    size_t curEpoch, size_t otherEpoch);
	static bool _dontAugmentDataSeries(const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                   const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                   const std::vector<std::vector<TypeTime>> &t, size_t l, size_t curEpoch,
	                                   size_t otherEpoch, TypeTime EpochBorder);
	static bool _epochIsIdentifiable(size_t L, const std::vector<TypeTime> &t,
	                                 const std::vector<size_t> &dataSeriesIndex,
	                                 const std::vector<bool> &AllDataPointsAreZero);
	static bool _epochIsNotIdentifiableBecauseOfZeroCounts(size_t L, const std::vector<TypeTime> &t,
	                                                       const std::vector<size_t> &dataSeriesIndex,
	                                                       const std::vector<size_t> &numDataPointsWithZeroCounts);

public:
	TOLSGamma(size_t NumEpochs, const std::vector<TypeTime> &EpochStarts);
	std::vector<double> estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints,
	                                         size_t SpeciesID);
};

//------------------------------------------
// TGamma
//------------------------------------------

class TGamma {
	// wraps parameters gamma and timesOfChange
protected:
	// stattools parameters
	TypeParamGamma *_gamma       = nullptr;
	TypeParamTOC *_timesOfChange = nullptr;

	// number of epochs
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<double> _rho;

	// Jeffrey's prior
	bool _noJP    = false;
	double _tryJP = 0.0;
	double _curJP = 0.0;

	// keep track of posterior summaries of gamma's
	std::vector<std::vector<size_t>> _summaryGammaPosterior;
	size_t _counter = 0;

	// initialization
	void _adjustTimesOfChange(const TUniqueContainer<TypeTime> &Timepoints);
	void _fillRho(const TUniqueContainer<TypeTime> &Timepoints);

public:
	TGamma(TypeParamGamma *Gamma, TypeParamTOC *TimesOfChange, const TUniqueContainer<TypeTime> &Timepoints);
	~TGamma() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t SpeciesID);
	void initializeJP(const TData &Data, const TModelBase &Model);

	std::vector<TypeTime> getEpochStarts(const TUniqueContainer<TypeTime> &Timepoints) const;

	double getLogPriorRatio(const TData &Data, const TModelBase &Model);

	bool updateSpecificIndex(size_t m);
	bool acceptOrReject(double logH, size_t m);

	void updateSummaryGammaPosterior();
	void resetSummaryGammaPosterior();
	void writeSummaryGammaPosterior(const std::string &Prefix);

	void setSimulatedGammaToZero();

	size_t numEpochs() const;
	size_t numTimepoints() const;
	TypeGamma gamma(size_t m) const;
	double rho(size_t row, size_t col) const;
	std::vector<TypePhi> phi() const;
	bool gammaIsUpdated() const;
};

//------------------------------------------
// TStochasticPrior
//------------------------------------------

class TStochasticPrior {
	// yellow box
protected:
	TypeParamLogPhi *_logPhi     = nullptr;
	TypeParamLogSigma *_logSigma = nullptr;

	// TGamma
	TGamma _gamma;

	// rho
	coretools::TMatrix<double> _stochasticRho;
	std::vector<double> _deltaTimeVec;

	// temporary variable: sum of rho * gamma
	std::vector<double> _newSumRhoGamma;
	std::vector<double> _oldSumRhoGamma;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<double>> _tryLogTransitionProbs;
	std::vector<std::vector<double>> _curLogTransitionProbs;

	// store determinitic phi
	std::vector<TypePhi> _newDeterministicPhi;
	std::vector<TypePhi> _oldDeterministicPhi;

	// fill rho and delta
	void _fillStochasticRho(const TUniqueContainer<TypeTime> &Timepoints);
	void _initializeSumRhoGamma();
	void _initializeTransitionProbabilities();

	// transition probabilities
	double _calculateLogTransitionProbability(size_t j, size_t k) const;
	coretools::probdist::TNormalDistr _getMeanSdTransitionProbability(size_t j, size_t k) const;

public:
	TStochasticPrior(TypeParamLogPhi *LogPhi, TypeParamLogSigma *LogSigma, TypeParamGamma *Gamma,
	                 TypeParamTOC *TimesOfChange, const TUniqueContainer<TypeTime> &Timepoints);
	~TStochasticPrior() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t SpeciesID);
	void initializeTempVariables(const TData &Data, const TModelBase &Model, bool CalculateJP);

	// update gamma
	void updateGamma(const TData &Data, const TModelBase &Model);
	double calculateLL() const;
	double calculateLLRatioForUpdateGamma();
	void swapTryCur();
	void swapTryCur(size_t j, size_t k);
	double calculateLogPriorRatio(size_t j, size_t k);

	// joint update of gamma
	std::vector<double> updateGamma_forLogPhiJointUpdate(size_t m);
	bool evaluateGamma_forLogPhiJointUpdate(size_t m, double LLRatio, const TData &Data, const TModelBase &Model);

	// write posterior summaries of gamma's
	void updateSummaryGammaPosterior();
	void resetSummaryGammaPosterior();
	void writeSummaryGammaPosterior(const std::string &Prefix);

	// simulations
	void setSimulatedGammaToZero();
	double sampleFromTransitionProbability(size_t j, size_t k) const;

	// getters
	TypePhi deterministicPhi(size_t k) const;
	std::vector<TypePhi> calculateDeterministicPhi() const;
	size_t numEpochs() const;
	bool gammaIsUpdated() const;
};

#endif // BIRP_TGAMMA_H
