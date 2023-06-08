//
// Created by caduffm on 4/13/23.
//

#ifndef BIRP_TMODELS_H
#define BIRP_TMODELS_H

#include "BirpTypes.h"
#include "TBirpConfiguration.h"
#include "TData.h"
#include "TDistributions.h"
#include "TGamma.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

class TGamma;

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
	static arma::mat _calcFisherInfo(const TGamma &Gamma, const std::vector<double> &p,
	                                 const std::vector<double> &sigma, TypeNu nu, const TLocations &loc);

protected:
	size_t _speciesIDinUniqueContainer = 0;

	void _initializeDistribution(const TData &Data, const TUniqueContainer<std::string> &Locations);
	void _estimateInitialParameters(const TData &Data);
	void _updateDistributionParameters(const TData &Data);
	void _simulateDistributionParameters(TData &Data, const std::vector<double> &AvgNumObsPerTimepoint,
	                                     bool SimulateCovariates);
	void _initializeLL(const TData &Data);

public:
	TModelBase(size_t SpeciesIDinUniqueContainer, const TData &Data,
	           const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn);
	virtual ~TModelBase() = default;

	virtual void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                        const TUniqueContainer<std::size_t> &TimePoints)                                      = 0;
	virtual void estimateInitial(const TData &Data, const TUniqueContainer<size_t> &Timepoints)                   = 0;
	virtual void initializeTempVariables(const TData &Data)                                                       = 0;
	virtual void update(const TData &Data)                                                                        = 0;
	virtual void simulate(TData &Data, const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) = 0;

	// calculate sigma
	std::vector<double> calculateSigma(const TLocations &loc) const;
	std::vector<double> calculateSigmaFromDeterministicPhi(const TLocations &loc) const;

	// abstract function for phi
	virtual TypePhi phi(size_t location_id, size_t timepoint_k) const = 0;
	virtual TypePhi deterministicPhi(size_t timepoint_k) const        = 0;

	double calculateJeffreyPrior(const TData &Data, const TGamma &Gamma) const;
	double calculateLLRatio(const TData &Data);
	double calculateLLRatio_perLocation(size_t j, const TData &Data);

	double calculateLLRatio_perMethod(size_t i, const TData &Data);

	void swapTryCur();
	void swapTryCur_perLocation(size_t j, const TData &Data);
	void swapTryCur_perMethod(size_t i, const TData &Data);

	double curLL(size_t ij) const;
	double sumCurLL() const;
	void setCurLL(size_t ij, double Value);
	void setCurLL(const std::vector<double> &Values);

	// write posterior summaries of gamma's
	virtual void resetSummaryGammaPosterior()                          = 0;
	virtual void writeSummaryGammaPosterior(const std::string &Prefix) = 0;

	virtual stattools::TUpdateBase *getPtrToJointUpdater() { return nullptr; }
};

//------------------------------------------
// TDeterministic
//------------------------------------------

class TDeterministic : public TModelBase {
protected:
	// gamma and times of change
	TGamma _gamma;

	// store temporary phi
	std::vector<TypePhi> _newPhi;
	std::vector<TypePhi> _oldPhi;

	void _fillPhi();

public:
	TDeterministic(size_t SpeciesIDinUniqueContainer, const TData &Data,
	               const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
	               stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	               stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TDeterministic() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                const TUniqueContainer<std::size_t> &TimePoints) override;
	void estimateInitial(const TData &Data, const TUniqueContainer<size_t> &Timepoints) override;
	void initializeTempVariables(const TData &Data) override;
	void update(const TData &Data) override;

	void simulate(TData &Data, const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) override;

	TypePhi phi(size_t location_id, size_t timepoint_id) const override;
	TypePhi deterministicPhi(size_t timepoint_k) const override;

	// write posterior summaries of gamma's
	void resetSummaryGammaPosterior() override;
	void writeSummaryGammaPosterior(const std::string &Prefix) override;
};

//------------------------------------------
// TStochastic
//------------------------------------------

class TStochastic : public TModelBase {
protected:
	// ptr to prior above
	TStochasticPrior _prior;

	// parameters: non-owning
	stattools::TParameterTyped<TypeLogPhi, 2> *_logPhi     = nullptr;
	stattools::TParameterTyped<TypeLogSigma, 1> *_logSigma = nullptr;

	std::unique_ptr<stattools::TUpdateUniqueIrregular<TypeLogPhi>> _jointUpdater;

	// initialization
	void _setPhiFromDeterministic();

	// update single logPhi
	void _updateLogPhi(size_t j, size_t k, const TData &Data);
	void _updateLogPhi(const TData &Data);

	// joint update logPhi and gamma
	void _jointUpdateLogPhiGamma(const TData &Data);

	// update logSigma
	void _updateLogSigma();

public:
	TStochastic(size_t SpeciesIDinUniqueContainer, const TData &Data, size_t NumLocations, size_t NumTimePoint,
	            const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
	            stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
	            stattools::TParameterTyped<TypeLogPhi, 2> *LogPhi,
	            stattools::TParameterTyped<TypeLogSigma, 1> *Logigma);
	~TStochastic() override = default;

	void initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
	                const TUniqueContainer<std::size_t> &TimePoints) override;
	void estimateInitial(const TData &Data, const TUniqueContainer<size_t> &Timepoints) override;
	void initializeTempVariables(const TData &Data) override;
	void update(const TData &Data) override;

	void simulate(TData &Data, const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) override;

	TypePhi phi(size_t location_id, size_t timepoint_k) const override;
	TypePhi deterministicPhi(size_t timepoint_k) const override;

	// write posterior summaries of gamma's
	void resetSummaryGammaPosterior() override;
	void writeSummaryGammaPosterior(const std::string &Prefix) override;

	stattools::TUpdateBase *getPtrToJointUpdater() override { return _jointUpdater.get(); }
};

#endif // BIRP_TMODELS_H
