#ifndef BIRP_TMODELS_H
#define BIRP_TMODELS_H

#include <vector>

#include "coretools/Math/TMatrix.h"
#include <armadillo>

#include "BirpDAG.h"
#include "BirpTypes.h"
#include "TBirpConfiguration.h"
#include "TData.h"
#include "TDistributions.h"
#include "TGamma.h"

class TPhi;

//---------------------------------------------------
// TModelBase
// pure abstract base class of TStochastic and TDeterministic
//---------------------------------------------------

class TModelBase {
public:
private:
	// temporary values of likelihood
	std::vector<double> _tryLL;
	std::vector<double> _curLL;

	// distribution: Poisson or negative binomial
	std::unique_ptr<TDistributionBase> _distr;

	// Fisher information
	arma::mat _calcRPR(size_t CI_ix, const TPhi &phi, double sum_sigma, const std::vector<double> &sigma, TypeNu nu,
	                   const TLocations &loc) const;

protected:
	size_t _speciesIDinUniqueContainer = 0;

	void _initializeDistribution(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                             TBirpPrior *BirpPrior);
	void _estimateInitialParameters(const TData &Data);
	void _updateDistributionParameters(const TData &Data);
	void _simulateDistributionParameters(TData &Data, bool use_n_bar, double n);
	void _simulateConditionedDistributionParameters(TData &Data);
	void _initializeLL(const TData &Data);
	void _updateTryLL(const TData &Data, size_t index_in_counts, size_t ix_method, size_t j_in_method);

public:
	TModelBase(size_t SpeciesIDinUniqueContainer, const TData &Data,
	           const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn);
	virtual ~TModelBase() = default;

	virtual void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                        const TUniqueContainer<TypeTime> &TimePoints, TBirpPrior *BirpPrior)  = 0;
	virtual void estimateInitial(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints) = 0;
	virtual void initializeTempVariables(const TData &Data)                                       = 0;
	virtual void update(const TData &Data)                                                        = 0;
	virtual void simulate(TData &Data, bool use_n_bar, double n)                                  = 0;
	virtual void simulateConditioned(TData &Data)                                                 = 0;

	// calculate sigma
	[[nodiscard]] std::vector<double> calculateSigma(const TLocations &loc) const;
	[[nodiscard]] std::vector<double> calculateSigmaFromDeterministicPhi(const TLocations &loc) const;

	// abstract function for phi
	[[nodiscard]] virtual TypePhi phi(size_t location_id, size_t timepoint_k, size_t CI_index) const = 0;
	[[nodiscard]] virtual TypePhi deterministicPhi(size_t timepoint_k, size_t CI_index) const        = 0;

	arma::mat calculateRPR_JP(const TData &Data, const TPhi &Phi, size_t CI_index) const;
	[[nodiscard]] double calculateJeffreyPrior(const TData &Data, const TPhi &Phi, size_t CI_index) const;
	double calculateLLRatio(const TData &Data);
	double calculateLLRatio_perLocation(size_t j, const TData &Data);
	double calculateLLRatio_perCIGroup(size_t CI_index, const TData &Data);
	double calculateLLRatio_perMethod(size_t i, const TData &Data);

	void swapTryCur();
	void swapTryCur_perLocation(size_t j, const TData &Data);
	void swapTryCur_perMethod(size_t i, const TData &Data);
	void swapTryCur_perCIGroup(size_t CI_index, const TData &Data);

	[[nodiscard]] double curLL(size_t ij) const;
	double curLLPerMethod(size_t i, const TData &Data);
	[[nodiscard]] double sumCurLL() const;
	void setCurLL(size_t ij, double Value);
	void setCurLLPerMethod(const std::vector<double> &Values, size_t i, const TData &Data);

	// write posterior summaries of gamma's and Delta's
	virtual void resetSummaryPosterior()                          = 0;
	virtual void writeSummaryPosterior(const std::string &Prefix) = 0;

	virtual stattools::TUpdateBase *getPtrToJointUpdater() { return nullptr; }
};

//------------------------------------------
// TDeterministic
//------------------------------------------

class TDeterministic : public TModelBase {
protected:
	// gamma and times of change
	TPhi _phi;

	// store temporary phi
	std::vector<std::vector<TypePhi>> _newPhi;
	std::vector<std::vector<TypePhi>> _oldPhi;

	void _fillPhi(size_t CI_ix);
	void _simulate();

	template<bool IsGamma> void _update(const TData &Data) {
		for (size_t ix = 0; ix < _phi.size<IsGamma>(); ix++) { // loop over all gamma / delta
			if (!_phi.updateSpecificIndex<IsGamma>(ix)) { continue; }

			double LLRatio = 0.0;
			for (const auto CI_index : _phi.getCIGroups<IsGamma>(ix)) {
				_oldPhi[CI_index] = _newPhi[CI_index];
				_fillPhi(CI_index);
				LLRatio += calculateLLRatio_perCIGroup(CI_index, Data);
			}
			const double logPriorRatio = _phi.getLogPriorRatio<IsGamma>(ix, Data, *this); // Jeffrey prior
			const double logH          = LLRatio + logPriorRatio;

			if (_phi.acceptOrReject<IsGamma>(logH, ix)) {
				// accepted
				for (const auto CI_index : _phi.getCIGroups<IsGamma>(ix)) { swapTryCur_perCIGroup(CI_index, Data); }
			} else {
				// reject
				for (const auto CI_index : _phi.getCIGroups<IsGamma>(ix)) { _newPhi[CI_index] = _oldPhi[CI_index]; }
			}
		}

		// update gamma and Delta summaries
		_phi.updateSummaryPosterior<IsGamma>();
	}

public:
	TDeterministic(size_t SpeciesIDinUniqueContainer, const TData &Data,
	               const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn, TypeParamGamma *Gamma,
	               TypeParamDelta *Delta, const TUniqueContainer<TypeTime> &Timepoints,
	               const TUniqueContainer<std::string> &CIGroupNames, std::string_view Prefix);
	~TDeterministic() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                const TUniqueContainer<TypeTime> &TimePoints, TBirpPrior *BirpPrior) override;
	void estimateInitial(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints) override;
	void initializeTempVariables(const TData &Data) override;
	void update(const TData &Data) override;

	void simulate(TData &Data, bool use_n_bar, double n) override;
	void simulateConditioned(TData &Data) override;

	[[nodiscard]] TypePhi phi(size_t location_id, size_t timepoint_id, size_t CI_index) const override;
	[[nodiscard]] TypePhi deterministicPhi(size_t timepoint_k, size_t CI_index) const override;

	// write posterior summaries of gamma's
	void resetSummaryPosterior() override;
	void writeSummaryPosterior(const std::string &Prefix) override;
};

//------------------------------------------
// TStochastic
//------------------------------------------

class TStochastic : public TModelBase {
protected:
	// ptr to prior above
	TStochasticPrior _prior;

	// parameters: non-owning
	TypeParamLogPhi *_logPhi     = nullptr;
	TypeParamLogSigma *_logSigma = nullptr;
	TypeParamLogZeta *_logZeta   = nullptr;

	std::unique_ptr<stattools::TUpdateUnique<TypeLogPhi, false>> _jointUpdater;

	// initialization
	void _setPhiFromDeterministic(const TData &Data);

	// update single logPhi
	void _updateLogPhi(size_t j, size_t k, const TData &Data);
	void _updateLogPhi(const TData &Data);

	// joint update logPhi and gamma
	void _jointUpdateLogPhiGamma(const TData &Data);

	// joint update of logPhi
	void _shiftLogPhis(size_t j, size_t k, double diff);
	void _jointUpdateLogPhi(size_t j, size_t k, const TData &Data);
	void _jointUpdateLogPhi(const TData &Data);

	// update logSigma and logZeta
	void _updateLogSigma(const TData &Data);
	void _updateLogZeta(size_t i, const TData &Data);

	// simulate
	void _simulate(TData &Data);

public:
	TStochastic(size_t SpeciesIDinUniqueContainer, const TData &Data, size_t NumLocations, size_t NumTimePoint,
	            const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn, TypeParamGamma *Gamma, TypeParamDelta *Delta,
	            TypeParamLogPhi *LogPhi, TypeParamLogSigma *LogSigma, TypeParamLogZeta *LogZeta,
	            const TUniqueContainer<TypeTime> &Timepoints, const TUniqueContainer<std::string> &CIGroupNames,
	            std::string_view Prefix);
	~TStochastic() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                const TUniqueContainer<TypeTime> &TimePoints, TBirpPrior *BirpPrior) override;
	void estimateInitial(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints) override;
	void initializeTempVariables(const TData &Data) override;
	void update(const TData &Data) override;

	void simulate(TData &Data, bool use_n_bar, double n) override;
	void simulateConditioned(TData &Data) override;

	[[nodiscard]] TypePhi phi(size_t location_id, size_t timepoint_k, size_t CI_index) const override;
	[[nodiscard]] TypePhi deterministicPhi(size_t timepoint_k, size_t CI_index) const override;

	// write posterior summaries of gamma's
	void resetSummaryPosterior() override;
	void writeSummaryPosterior(const std::string &Prefix) override;

	stattools::TUpdateBase *getPtrToJointUpdater() override { return _jointUpdater.get(); }
};

#endif // BIRP_TMODELS_H
