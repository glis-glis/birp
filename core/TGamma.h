#ifndef BIRP_TGAMMA_H
#define BIRP_TGAMMA_H

#include "coretools/Math/TMatrix.h"
#include <armadillo>
#include <vector>

#include "BirpDAG.h"
#include "BirpTypes.h"
#include "TUniqueContainer.h"

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
	size_t _numGamma;

	// BACI
	std::vector<std::vector<size_t>> _gamma_ix_per_CIgroup_and_epoch;
	std::vector<std::vector<size_t>> _CI_indices_per_gamma;

	// number of epochs
	size_t _numEpochs = 0;
	std::vector<TypeTime> _epochStarts;

	bool _epochIsIdentifiable(size_t curGamma, const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                          const std::vector<std::vector<size_t>> &numDataPointsWithZeroCounts) const;

	bool _timePointIsRelevant(const TUniqueContainer<TypeTime> &Timepoints, size_t t, size_t e) const;
	auto _loopOverTimepoints(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t SpeciesID,
	                         size_t i, size_t j_in_method, size_t e, size_t L, std::vector<double> &y_ix,
	                         std::vector<TypeTime> &t_ix, std::vector<size_t> &dataSeriesIndex) const;
	auto _getNextGammaToInitialize(const std::vector<bool> &GammaInitialized,
	                               const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                               const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                               const std::vector<std::vector<size_t>> &numDataPointsWithZeroCounts) const;
	auto _get_Y_ix(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t ix, size_t SpeciesID) const;
	static arma::vec _doOLS(const std::vector<double> &Y_m, const std::vector<TypeTime> &T_m,
	                        const std::vector<size_t> &DataSeriesIndex, size_t L);
	void _augmentYWithOLSEstimatesStart(size_t curGamma, size_t CI_index, size_t epoch, const arma::vec &Beta,
	                                    std::vector<std::vector<double>> &y, std::vector<std::vector<TypeTime>> &t,
	                                    std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                    bool augment_zero) const;
	void _augmentYWithOLSEstimatesStart_All(size_t curGamma, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                        std::vector<std::vector<TypeTime>> &t,
	                                        std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                        bool augment_zero) const;
	void _augmentYWithOLSEstimatesEnd(size_t curGamma, size_t CI_index, size_t epoch, const arma::vec &Beta,
	                                  std::vector<std::vector<double>> &y, std::vector<std::vector<TypeTime>> &t,
	                                  std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                  bool augment_zero) const;
	void _augmentYWithOLSEstimatesEnd_All(size_t curGamma, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                      std::vector<std::vector<TypeTime>> &t,
	                                      std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                      std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                      std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                      bool augment_zero) const;
	static void
	_augmentDataSeriesIndexWithOLSEstimates(std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                        size_t l, std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                        size_t curGamma, size_t otherGamma);
	static bool _dontAugmentDataSeries(const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                   const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                   const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                   const std::vector<std::vector<TypeTime>> &t, size_t CI_index, size_t epoch,
	                                   size_t l, size_t curGamma, size_t otherGamma, TypeTime EpochBorder);
	static auto _getScore(size_t l, const std::vector<size_t> &dataSeriesIndex,
	                      const std::vector<size_t> &numDataPointsWithZeroCounts);
	static size_t _getNumNonZeroDataPointsInRange(size_t l, const std::vector<size_t> &dataSeriesIndex);

public:
	TOLSGamma(size_t NumGamma, size_t NumEpochs, const std::vector<TypeTime> &EpochStarts,
	          const std::vector<std::vector<size_t>> &Gamma_ix_per_CIgroup_and_epoch,
	          const std::vector<std::vector<size_t>> &CI_indices_per_gamma);
	std::vector<double> estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints,
	                                         size_t SpeciesID, const TUniqueContainer<std::string> &GammaNames);
};

//------------------------------------------
// TPhi
//------------------------------------------

class TPhi {
	// wraps parameters gamma and timesOfChange
protected:
	// stattools parameters
	TypeParamGamma *_gamma = nullptr;
	TypeParamDelta *_delta = nullptr;

	// species id
	size_t _species_id = 0;

	// times of change
	std::vector<TypeTime> _timesOfChange;
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<double> _rho;
	coretools::TMatrix<double> _ind_delta;

	// Jeffrey's prior
	bool _noJP = false;
	std::vector<double> _tryJP;
	std::vector<double> _curJP;

	// BACI
	size_t _numCIGroups;
	std::vector<std::vector<size_t>> _gamma_ix_per_CIgroup_and_epoch;
	TUniqueContainer<std::string> _gamma_names;
	std::vector<std::vector<size_t>> _CI_indices_per_gamma;
	std::vector<std::vector<size_t>> _epoch_indices_per_gamma;

	std::vector<std::vector<size_t>> _delta_ix_per_CIgroup_and_epoch;
	TUniqueContainer<std::string> _delta_names;
	std::vector<std::vector<size_t>> _CI_indices_per_delta;

	// keep track of posterior summaries of gamma's
	std::vector<std::vector<size_t>> _summaryGammaPosterior;
	size_t _counter = 0;

	// initialization
	void _fillTimesOfChange(const TUniqueContainer<TypeTime> &Timepoints, std::string_view Prefix);
	void _fillRho(const TUniqueContainer<TypeTime> &Timepoints);
	void _fillIndicatorDelta(const TUniqueContainer<TypeTime> &Timepoints);
	void _readGrouping(bool is_gamma, std::string_view Name, const TUniqueContainer<TypeTime> &Timepoints,
	                   const TUniqueContainer<std::string> &CIGroupNames, size_t &numCIGroups,
	                   std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch, TUniqueContainer<std::string> &names,
	                   std::vector<std::vector<size_t>> &CI_indices, std::string_view Prefix) const;
	void _fillIndicatorsBA(bool is_gamma, std::string_view Name, const TUniqueContainer<TypeTime> &Timepoints,
	                       std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                       TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices) const;
	void _fillIndicatorsCI(std::string_view Name, const TUniqueContainer<std::string> &CIGroupNames,
	                       std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                       TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices) const;
	void _fillIndicatorsBACI(std::string_view Filename, std::string_view Name,
	                         const TUniqueContainer<std::string> &CIGroupNames,
	                         std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                         TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices) const;

public:
	TPhi(TypeParamGamma *Gamma, TypeParamDelta *Delta, const TUniqueContainer<TypeTime> &Timepoints,
	     const TUniqueContainer<std::string> &CIGroupNames, size_t SpeciesID, std::string_view Prefix);
	~TPhi() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints);
	void initializeJP(const TData &Data, const TModelBase &Model);

	std::vector<TypeTime> getEpochStarts(const TUniqueContainer<TypeTime> &Timepoints) const;

	double getLogPriorRatio(const TData &Data, const TModelBase &Model, size_t CI_ix);

	template<bool IsGamma> bool updateSpecificIndex(size_t ix) {
		if constexpr (IsGamma) {
			if (_gamma->isUpdated()) {
				_gamma->propose(coretools::TRange(ix));
				return true;
			}
			return false;
		} else {
			if (_delta->isUpdated()) {
				_delta->propose(coretools::TRange(ix));
				return true;
			}
			return false;
		}
	}

	template<bool IsGamma> bool acceptOrReject(double logH, size_t ix) {
		if constexpr (IsGamma) {
			if (_gamma->acceptOrReject(logH, coretools::TRange(ix))) {
				for (const auto CI_ix : _CI_indices_per_gamma[ix]) { _curJP[CI_ix] = _tryJP[CI_ix]; }
				return true;
			}
			return false;
		} else {
			if (_delta->acceptOrReject(logH, coretools::TRange(ix))) {
				for (const auto CI_ix : _CI_indices_per_delta[ix]) { _curJP[CI_ix] = _tryJP[CI_ix]; }
				return true;
			}
			return false;
		}
	}

	void updateSummaryGammaPosterior();
	void resetSummaryGammaPosterior();
	void writeSummaryGammaPosterior(const std::string &Prefix);

	void setSimulatedGammaToZero();

	size_t numEpochs() const;
	size_t numCIGroups() const;
	size_t numTimepoints() const;
	size_t numGamma() const;
	size_t numDelta() const;
	size_t species_id() const;
	template<bool IsGamma> size_t size() const {
		if constexpr (IsGamma) { return numGamma(); }
		return numDelta();
	}
	template<bool IsGamma> const std::vector<size_t> &getCIGroups(size_t gamma_ix) const {
		if constexpr (IsGamma) { return _CI_indices_per_gamma[gamma_ix]; }
		return _CI_indices_per_delta[gamma_ix];
	}

	TypeGamma gamma(size_t m, size_t CI_ix) const;
	TypeDelta delta(size_t epoch, size_t CI_ix) const;
	double rho(size_t row, size_t col) const;
	double indDelta(size_t row, size_t col) const;
	std::vector<TypePhi> phi(size_t CI_ix) const;
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
	std::vector<std::vector<double>> _newSumRhoGamma;
	std::vector<std::vector<double>> _oldSumRhoGamma;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<double>> _tryLogTransitionProbs;
	std::vector<std::vector<double>> _curLogTransitionProbs;

	// store deterministic phi
	std::vector<std::vector<TypePhi>> _newDeterministicPhi;
	std::vector<std::vector<TypePhi>> _oldDeterministicPhi;

	// fill rho and delta
	void _fillStochasticRho(const TUniqueContainer<TypeTime> &Timepoints);
	void _initializeSumRhoGamma(size_t CI_index);
	void _initializeTransitionProbabilities(const TData &Data);

	// transition probabilities
	double _calculateLogTransitionProbability(size_t j, size_t k, size_t CI_index) const;
	coretools::probdist::TNormalDistr _getMeanSdTransitionProbability(size_t j, size_t k, size_t CI_index) const;

public:
	TStochasticPrior(TypeParamLogPhi *LogPhi, TypeParamLogSigma *LogSigma, TypeParamGamma *Gamma, TypeParamDelta *Delta,
	                 const TUniqueContainer<TypeTime> &Timepoints,
	                 const TUniqueContainer<std::string> &CIGroupNamesGamma, size_t SpeciesID, std::string_view Prefix);
	~TStochasticPrior() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints);
	void initializeTempVariables(const TData &Data, const TModelBase &Model, bool CalculateJP);

	// update gamma
	void updateGamma(const TData &Data, const TModelBase &Model);
	double calculateLL(const TData &Data, size_t CI_index) const;
	double calculateLLRatioForUpdateGamma(const TData &Data, size_t CI_index);
	void swapTryCur();
	void swapTryCur(size_t ix, const TData &Data);
	void swapTryCur(size_t j, size_t k);
	double calculateLogPriorRatio(size_t j, size_t k, size_t CI_index);

	// joint update of gamma
	void updateGamma_forLogPhiJointUpdate(size_t ix);
	std::vector<double> getDelta_forLogPhiJointUpdate(size_t CI_index);
	bool evaluateGamma_forLogPhiJointUpdate(size_t m, double LLRatio, const TData &Data, const TModelBase &Model);

	// write posterior summaries of gamma's
	void updateSummaryGammaPosterior();
	void resetSummaryGammaPosterior();
	void writeSummaryGammaPosterior(const std::string &Prefix);

	// simulations
	void setSimulatedGammaToZero();
	double sampleFromTransitionProbability(size_t j, size_t k, size_t CI_index) const;

	// getters
	TypePhi deterministicPhi(size_t k, size_t CI_index) const;
	std::vector<TypePhi> calculateDeterministicPhi(size_t CI_index) const;
	size_t numEpochs() const;
	bool gammaIsUpdated() const;
	size_t numCIGroups() const;
	size_t size() const;
	const std::vector<size_t> &getCIGroups(size_t gamma_ix) const;
};

#endif // BIRP_TGAMMA_H
