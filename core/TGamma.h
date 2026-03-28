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

	bool _timePointIsRelevant(const TUniqueContainer<TypeTime> &Timepoints, size_t t, size_t e) const;
	bool _timePointIsRelevantDelta(const TUniqueContainer<TypeTime> &Timepoints, size_t t, size_t e) const;
	auto _loopOverTimepoints(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t SpeciesID,
	                         size_t i, size_t j_in_method, size_t e, size_t CI_index, size_t L,
	                         std::vector<double> &y_ix, std::vector<TypeTime> &t_ix,
	                         std::vector<size_t> &dataSeriesIndex,
	                         const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                         const std::vector<bool> &update_delta) const;
	auto _getNextGammaToInitialize(const std::vector<bool> &GammaInitialized,
	                               const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                               const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                               const std::vector<std::vector<size_t>> &numDataPointsWithZeroCounts,
	                               const std::vector<bool> &update_gamma) const;
	auto _get_Y_ix(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t ix, size_t SpeciesID,
	               const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	               const std::vector<bool> &update_delta) const;
	static arma::vec _doOLS(const std::vector<double> &Y_m, const std::vector<TypeTime> &T_m,
	                        const std::vector<size_t> &DataSeriesIndex, size_t L);
	void _augmentYWithOLSEstimatesStart(size_t curGamma, size_t CI_index, size_t epoch, const arma::vec &Beta,
	                                    std::vector<std::vector<double>> &y, std::vector<std::vector<TypeTime>> &t,
	                                    std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                    std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                    size_t epoch_ix_delta,
	                                    const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                    const std::vector<bool> &update_delta, bool augment_zero) const;
	void _augmentYWithOLSEstimatesStart_All(size_t curGamma, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                        std::vector<std::vector<TypeTime>> &t,
	                                        std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                        const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                        const std::vector<bool> &update_delta, bool augment_zero) const;
	void _augmentYWithOLSEstimatesEnd(size_t curGamma, size_t CI_index, size_t epoch, const arma::vec &Beta,
	                                  std::vector<std::vector<double>> &y, std::vector<std::vector<TypeTime>> &t,
	                                  std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                  std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                  size_t epoch_ix_delta,
	                                  const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                  const std::vector<bool> &update_delta, bool augment_zero) const;
	void _augmentYWithOLSEstimatesEnd_All(size_t curGamma, const arma::vec &Beta, std::vector<std::vector<double>> &y,
	                                      std::vector<std::vector<TypeTime>> &t,
	                                      std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                      std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                      std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                      const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                      const std::vector<bool> &update_delta, bool augment_zero) const;
	static void
	_augmentDataSeriesIndexWithOLSEstimates(std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                        std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                        size_t l, std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                        size_t curGamma, size_t otherGamma);
	static bool _dontAugmentDataSeries(const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
	                                   const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_CI_epoch,
	                                   const std::vector<std::vector<size_t>> &dataSeriesIndex,
	                                   const std::vector<std::vector<TypeTime>> &t, size_t CI_index, size_t epoch,
	                                   size_t l, size_t curGamma, size_t otherGamma, TypeTime EpochBorder,
	                                   size_t epoch_ix_delta,
	                                   const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                   const std::vector<bool> &update_delta);
	static auto _getScore(size_t l, const std::vector<size_t> &dataSeriesIndex,
	                      const std::vector<size_t> &numDataPointsWithZeroCounts);
	static size_t _getNumNonZeroDataPointsInRange(size_t l, const std::vector<size_t> &dataSeriesIndex);
	void _checkIfDeltaIsIdentifiable(const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                 const std::vector<bool> &update_delta, const TData &Data,
	                                 const TUniqueContainer<TypeTime> &Timepoints) const;

public:
	TOLSGamma(size_t NumGamma, size_t NumEpochs, const std::vector<TypeTime> &EpochStarts,
	          const std::vector<std::vector<size_t>> &Gamma_ix_per_CIgroup_and_epoch,
	          const std::vector<std::vector<size_t>> &CI_indices_per_gamma);
	std::vector<double> estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints,
	                                         size_t SpeciesID, const TUniqueContainer<std::string> &GammaNames,
	                                         const std::vector<bool> &update_gamma,
	                                         const std::vector<std::vector<size_t>> &delta_ix_per_CIgroup_and_epoch,
	                                         const std::vector<bool> &update_delta);
};

//------------------------------------------
// TPhi
//------------------------------------------
enum class Changes { only_gamma, only_Delta, both };

class TPhi {
	// wraps parameters gamma and timesOfChange
protected:
	// stattools parameters
	TypeParamGamma *_gamma = nullptr;
	TypeParamDelta *_delta = nullptr;

	// which parameters to infer
	Changes _change;

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
	std::vector<arma::mat> _newRPR;
	std::vector<arma::mat> _oldRPR;
	double _newJP;
	double _oldJP;

	// BACI
	size_t _numCIGroups;

	std::vector<std::vector<size_t>> _gamma_ix_per_CIgroup_and_epoch;
	TUniqueContainer<std::string> _gamma_names;
	std::vector<std::vector<size_t>> _CI_indices_per_gamma;
	std::vector<bool> _update_gamma;

	std::vector<std::vector<size_t>> _delta_ix_per_CIgroup_and_epoch;
	TUniqueContainer<std::string> _delta_names;
	std::vector<std::vector<size_t>> _CI_indices_per_delta;
	std::vector<bool> _update_delta;

	// keep track of posterior summaries of gamma's and Deltas
	std::vector<std::vector<size_t>> _summaryGammaPosterior;
	std::vector<std::vector<size_t>> _summaryDeltaPosterior;
	size_t _counterGamma = 0;
	size_t _counterDelta = 0;

	// initialization
	void _read_change_to_infer();
	void _read_BACI_change(const TUniqueContainer<TypeTime> &Timepoints,
	                       const TUniqueContainer<std::string> &CIGroupNames, std::string_view Prefix);
	void _fillTimesOfChange(const TUniqueContainer<TypeTime> &Timepoints);
	void _fillRho(const TUniqueContainer<TypeTime> &Timepoints);
	void _fillIndicatorDelta(const TUniqueContainer<TypeTime> &Timepoints);
	void _readGrouping(bool is_gamma, std::string_view ParamName, std::string_view ArgNameDesign,
	                   const TUniqueContainer<TypeTime> &Timepoints, const TUniqueContainer<std::string> &CIGroupNames,
	                   std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch, TUniqueContainer<std::string> &names,
	                   std::vector<std::vector<size_t>> &CI_indices, std::vector<bool> &update,
	                   std::string_view Prefix) const;
	void _fillIndicatorsBA(bool is_gamma, std::string_view Name, const TUniqueContainer<TypeTime> &Timepoints,
	                       std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                       TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices) const;
	void _fillIndicatorsCI(std::string_view Name, const TUniqueContainer<std::string> &CIGroupNames,
	                       std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                       TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices) const;
	void _fillIndicatorsBACI(std::string_view Filename, const TUniqueContainer<std::string> &CIGroupNames,
	                         std::vector<std::vector<size_t>> &ix_per_CIgroup_and_epoch,
	                         TUniqueContainer<std::string> &names, std::vector<std::vector<size_t>> &CI_indices, std::string_view Name) const;

	template<class T> bool _update(size_t ix, T *Param, const std::vector<bool> &allow_update) {
		if (Param->isUpdated() && allow_update[ix]) {
			Param->propose(coretools::TRange(ix));
			return true;
		}
		return false;
	}

	template<class T>
	bool _acceptOrReject(double logH, size_t ix, T *Param, const std::vector<std::vector<size_t>> &CI_indices,
	                     const std::vector<bool> &allow_update) {
		if (!allow_update[ix]) { return true; }

		if (Param->acceptOrReject(logH, coretools::TRange(ix))) {
			// accepted -> new RPR and newJP are already correct, no need to do something
			return true;
		}
		// rejected -> reset new = old
		for (const auto CI_ix : CI_indices[ix]) { _newRPR[CI_ix] = _oldRPR[CI_ix]; }
		_newJP = _oldJP;
		return false;
	}

	template<class T>
	size_t _updateSummaryPosterior(size_t counter, T *Param, std::vector<std::vector<size_t>> &summaryPosterior) const {
		for (size_t ix_1 = 0; ix_1 < Param->size(); ++ix_1) {     // rows
			for (size_t ix_2 = 0; ix_2 < Param->size(); ++ix_2) { // column
				if (ix_1 == ix_2) {
					// on the diagonal: count if Param_ix_1 > 0
					if (Param->value(ix_1) >= 0.0) { ++summaryPosterior[ix_1][ix_2]; }
				} else {
					// off the diagonal: count if Param_ix_1 > Param_ix_2
					if (Param->value(ix_1) > Param->value(ix_2)) { ++summaryPosterior[ix_1][ix_2]; }
				}
			}
		}
		return counter + 1;
	}

	template<class T>
	void _writeSummaryPosterior(const std::string &filename, const TUniqueContainer<std::string> &names, T *Param,
	                            const std::vector<std::vector<size_t>> &summaryPosterior, size_t counter) {
		auto header = names.vec();
		header.insert(header.begin(), "row_bigger_than_column");
		coretools::TOutputMaybeRcppFile file(filename, header);
		for (size_t ix_1 = 0; ix_1 < Param->size(); ++ix_1) {
			file << names[ix_1];
			for (size_t ix_2 = 0; ix_2 < Param->size(); ++ix_2) {
				file << (double)summaryPosterior[ix_1][ix_2] / (double)counter;
			}
			file.endln();
		}
	}

	double _getLogPriorRatio(const std::vector<size_t> &CI_indices, const TData &Data, const TModelBase &Model,
	                         bool allow_update);
	double _calculateJP() const;

public:
	TPhi(TypeParamGamma *Gamma, TypeParamDelta *Delta, const TUniqueContainer<TypeTime> &Timepoints,
	     const TUniqueContainer<std::string> &CIGroupNames, size_t SpeciesID, std::string_view Prefix);
	~TPhi() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialValues(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints);
	void initializeJP(const TData &Data, const TModelBase &Model);

	std::vector<TypeTime> getEpochStarts(const TUniqueContainer<TypeTime> &Timepoints) const;

	template<bool IsGamma> bool updateSpecificIndex(size_t ix) {
		if constexpr (IsGamma) {
			return _update(ix, _gamma, _update_gamma);
		} else {
			return _update(ix, _delta, _update_delta);
		}
	}

	template<bool IsGamma> bool acceptOrReject(double logH, size_t ix) {
		if constexpr (IsGamma) { return _acceptOrReject(logH, ix, _gamma, _CI_indices_per_gamma, _update_gamma); }
		return _acceptOrReject(logH, ix, _delta, _CI_indices_per_delta, _update_delta);
	}

	template<bool IsGamma> double getLogPriorRatio(size_t ix, const TData &Data, const TModelBase &Model) {
		// -> prior on gamma and Delta: Jeffrey prior
		if (_noJP) { return 0.0; }
		if constexpr (IsGamma) { return _getLogPriorRatio(_CI_indices_per_gamma[ix], Data, Model, _update_gamma[ix]); }
		return _getLogPriorRatio(_CI_indices_per_delta[ix], Data, Model, _update_delta[ix]);
	}

	template<bool IsGamma> void updateSummaryPosterior() {
		if constexpr (IsGamma) {
			_counterGamma = _updateSummaryPosterior(_counterGamma, _gamma, _summaryGammaPosterior);
		} else {
			_counterDelta = _updateSummaryPosterior(_counterDelta, _delta, _summaryDeltaPosterior);
		}
	}

	template<bool IsGamma> void resetSummaryPosterior() {
		if constexpr (IsGamma) {
			_counterGamma = 0;
			for (auto &m : _summaryGammaPosterior) { std::fill(m.begin(), m.end(), 0); }
		} else {
			_counterDelta = 0;
			for (auto &m : _summaryDeltaPosterior) { std::fill(m.begin(), m.end(), 0); }
		}
	}

	template<bool IsGamma> void writeSummaryPosterior(const std::string &Prefix) {
		if constexpr (IsGamma) {
			_writeSummaryPosterior(Prefix + "_gammaSummaries.txt", _gamma_names, _gamma, _summaryGammaPosterior,
			                       _counterGamma);
		} else {
			_writeSummaryPosterior(Prefix + "_DeltaSummaries.txt", _delta_names, _delta, _summaryDeltaPosterior,
			                       _counterDelta);
		}
	}

	void setSimulatedToZero();

	Changes change() const;
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
	template<bool IsGamma> const std::vector<size_t> &getCIGroups(size_t ix) const {
		if constexpr (IsGamma) { return _CI_indices_per_gamma[ix]; }
		return _CI_indices_per_delta[ix];
	}
	template<bool IsGamma> size_t getCIGroup(size_t CI_group_ix, size_t epoch) const {
		if constexpr (IsGamma) { return _CI_indices_per_gamma[CI_group_ix][epoch]; }
		return _CI_indices_per_delta[CI_group_ix][epoch];
	}
	template<bool IsGamma> size_t ix_per_CI_group_and_epoch(size_t CI_group_ix, size_t epoch) const {
		if constexpr (IsGamma) { return _gamma_ix_per_CIgroup_and_epoch[CI_group_ix][epoch]; }
		return _delta_ix_per_CIgroup_and_epoch[CI_group_ix][epoch];
	}

	TypeGamma gamma(size_t m, size_t CI_ix) const;
	TypeDelta delta(size_t epoch, size_t CI_ix) const;
	double rho(size_t row, size_t col) const;
	double indDelta(size_t row, size_t col) const;
	std::vector<TypePhi> phi(size_t CI_ix) const;
	bool gammaIsUpdated() const;
	bool deltaIsUpdated() const;
	const std::vector<std::string> &delta_names() const;
};

//------------------------------------------
// TStochasticPrior
//------------------------------------------

class TStochasticPrior {
	// yellow box
protected:
	TypeParamLogPhi *_logPhi     = nullptr;
	TypeParamLogSigma *_logSigma = nullptr;
	TypeParamLogZeta *_logZeta   = nullptr;

	// TGamma
	TPhi _phi;

	// rho
	coretools::TMatrix<double> _stochasticRho;
	coretools::TMatrix<double> _stochasticIndDelta;
	std::vector<double> _deltaTimeVec;

	// temporary variable: sum of rho * gamma + I * Delta
	std::vector<std::vector<double>> _newSumRhoGammaDelta;
	std::vector<std::vector<double>> _oldSumRhoGammaDelta;

	// temporary variable: sum of I * Zeta
	std::vector<std::vector<double>> _newSumZeta;
	std::vector<std::vector<double>> _oldSumZeta;

	// store log transition probabilities: per location and timepoint
	std::vector<std::vector<double>> _tryLogTransitionProbs;
	std::vector<std::vector<double>> _curLogTransitionProbs;

	// store deterministic phi
	std::vector<std::vector<TypePhi>> _newDeterministicPhi;
	std::vector<std::vector<TypePhi>> _oldDeterministicPhi;

	// fill rho and delta
	void _fillStochasticRho(const TUniqueContainer<TypeTime> &Timepoints);
	void _fillStochasticIndDelta(const TUniqueContainer<TypeTime> &Timepoints);
	void _initializeSumRhoGammaDelta(size_t CI_index);
	void _initializeSumZeta(size_t CI_index);
	void _initializeTransitionProbabilities(const TData &Data);

	// transition probabilities
	double _calculateLogTransitionProbability(size_t j, size_t k, size_t CI_index) const;
	coretools::probdist::TNormalDistr _getMeanSdTransitionProbability(size_t j, size_t k, size_t CI_index) const;
	void _swapTryCurGammaDelta(size_t CI_index, const TData &Data);

public:
	TStochasticPrior(TypeParamLogPhi *LogPhi, TypeParamLogSigma *LogSigma, TypeParamLogZeta *LogZeta,
	                 TypeParamGamma *Gamma, TypeParamDelta *Delta, const TUniqueContainer<TypeTime> &Timepoints,
	                 const TUniqueContainer<std::string> &CIGroupNamesGamma, size_t SpeciesID, std::string_view Prefix);
	~TStochasticPrior() = default;

	void initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior);
	void estimateInitialValues(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints);
	void initializeTempVariables(const TData &Data, const TModelBase &Model, bool CalculateJP);

	// update gamma/delta
	double calculateLL(const TData &Data, size_t CI_index) const;
	double calculateLLRatioForUpdateGammaDelta(const TData &Data, size_t CI_index);
	void swapTryCurSigma();
	void swapTryCurZeta(size_t ix, const TData &Data);
	void swapTryCurLogPhi(size_t j, size_t k);

	template<bool IsGamma> void updateGammaDelta(const TData &Data, const TModelBase &Model) {
		for (size_t ix = 0; ix < _phi.size<IsGamma>(); ix++) { // update gamma or delta one by one
			if (!_phi.updateSpecificIndex<IsGamma>(ix)) { continue; }

			double LLRatio = 0.0;
			for (auto CI_index : _phi.getCIGroups<IsGamma>(ix)) {
				// loop over all CI groups relevant for that gamma/delta
				// remember current sum rho
				_oldSumRhoGammaDelta[CI_index] = _newSumRhoGammaDelta[CI_index];
				_initializeSumRhoGammaDelta(CI_index);
				// remember current phi
				_oldDeterministicPhi[CI_index] = _newDeterministicPhi[CI_index];
				_newDeterministicPhi[CI_index] = _phi.phi(CI_index);
				LLRatio += calculateLLRatioForUpdateGammaDelta(Data, CI_index);
			}

			const double logPriorRatio = _phi.getLogPriorRatio<IsGamma>(ix, Data, Model); // Jeffrey prior
			const double logH          = LLRatio + logPriorRatio;

			if (_phi.acceptOrReject<IsGamma>(logH, ix)) { // already swaps JP (if necessary)
				// accepted
				swapTryCurGammaDelta<IsGamma>(ix, Data); // -> _curLogTransitionProbs = _tryLogTransitionProbs;
			} else {
				// rejected
				for (const auto CI_index : _phi.getCIGroups<IsGamma>(ix)) {
					_newSumRhoGammaDelta[CI_index] = _oldSumRhoGammaDelta[CI_index];
					_newDeterministicPhi[CI_index] = _oldDeterministicPhi[CI_index];
				}
			}
		}

		// update summaries of gamma and Delta posterior after update has finished
		updateSummaryPosterior<IsGamma>();
	}

	// update of zeta
	void updateSumZeta(size_t ix_zeta);
	void rejectSumZeta(size_t ix_zeta);

	template<bool IsGamma> void swapTryCurGammaDelta(size_t ix, const TData &Data) {
		// used when gamma or delta were updated
		for (const auto CI_index : _phi.getCIGroups<IsGamma>(ix)) { _swapTryCurGammaDelta(CI_index, Data); }
	}

	double calculateLogPriorRatio(size_t j, size_t k, size_t CI_index);

	// joint update of gamma
	template<bool IsGamma> bool update_forLogPhiJointUpdate(size_t ix) { return _phi.updateSpecificIndex<IsGamma>(ix); }
	std::vector<double> getDiff_forLogPhiJointUpdate(size_t CI_index);
	template<bool IsGamma>
	bool evaluate_forLogPhiJointUpdate(size_t ix, double LLRatio, const TData &Data, const TModelBase &Model) {
		double logH = LLRatio + _phi.getLogPriorRatio<IsGamma>(ix, Data, Model);

		if (_phi.acceptOrReject<IsGamma>(logH, ix)) { // already swaps JP (if necessary)
			// no need to swap transition probabilities: these didn't change!
			return true;
		}
		// rejected
		for (auto CI_index : _phi.getCIGroups<IsGamma>(ix)) {
			_newSumRhoGammaDelta[CI_index] = _oldSumRhoGammaDelta[CI_index];
			_newDeterministicPhi[CI_index] = _oldDeterministicPhi[CI_index];
		}

		return false;
	}

	// write posterior summaries of gamma's and Deltas
	template<bool IsGamma> void updateSummaryPosterior() { _phi.updateSummaryPosterior<IsGamma>(); }
	template<bool IsGamma> void resetSummaryGammaPosterior() { _phi.resetSummaryPosterior<IsGamma>(); }
	template<bool IsGamma> void writeSummaryGammaPosterior(const std::string &Prefix) {
		_phi.writeSummaryPosterior<IsGamma>(Prefix);
	}

	// simulations
	void setSimulatedToZero();
	double sampleFromTransitionProbability(size_t j, size_t k, size_t CI_index) const;

	// getters
	size_t numEpochs() const;
	size_t numCIGroups() const;

	const TPhi &phi() const;
	TypePhi deterministicPhi(size_t k, size_t CI_index) const;
	std::vector<TypePhi> calculateDeterministicPhi(size_t CI_index) const;

	Changes change() const;
	template<bool IsGamma> bool isUpdated() const {
		if constexpr (IsGamma) {
			return _phi.gammaIsUpdated();
		} else {
			return _phi.deltaIsUpdated();
		}
	}

	size_t numGamma() const;
	size_t numDelta() const;
	template<bool IsGamma> size_t size() const { return _phi.size<IsGamma>(); }
	template<bool IsGamma> const std::vector<size_t> &getCIGroups(size_t ix) const {
		return _phi.getCIGroups<IsGamma>(ix);
	}
	const std::vector<std::string> &delta_names() const;
	double zeta(size_t epoch, size_t CI_ix) const;
};

#endif // BIRP_TGAMMA_H
