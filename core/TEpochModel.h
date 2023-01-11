
//
// Created by Liam on 23/08/22.
//

#ifndef BIRP_TEPOCHMODEL_H
#define BIRP_TEPOCHMODEL_H

#include "BirpTypes.h"
#include "TSpecies.h"
#include "coretools/Math/TMatrix.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <armadillo>
#include <vector>

//--------------------------------
// TEpochModel
//--------------------------------

class TEpochModel {
protected:
	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1> *_gamma        = nullptr;
	stattools::TParameterTyped<TypeTime, 1> *_timesOfChange = nullptr;

	// number of epochs
	size_t _numEpochs = 0;

	// rho
	coretools::TMatrix<size_t> _rho;

	// temporary values, stored for speed
	// TODO: make sure old LL and old JP refer to correct species in case of multi species inference
	double _oldLL = 0.0;
	double _oldJP = 0.0;

	// likelihood
	[[nodiscard]] std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const;
	[[nodiscard]] double _calculateLogPhi_k(size_t k) const;
	[[nodiscard]] std::vector<double> _calculatePhi() const;
	arma::mat _calcFisherInfo(const std::vector<double> &p, const std::vector<double> &sigma, TypeCounts nu) const;
	std::pair<double, double> _calculateLLandJP(const TSpecies &Data) const;
	void _fillRho(const TUniqueContainer<size_t> &Timepoints);

	// simulation
	void _setSimulatedGammaToZero();

	virtual double _calculateLL(const std::vector<double> &sigma, const TDataSeries &lm, const std::vector<double> &p) const = 0;

public:
	TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	            stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
	~TEpochModel() = default;
	virtual void initializeInferred(const TUniqueContainer<size_t> &Timepoints,
	                                const TUniqueContainer<std::string> &Locations,
	                                const TUniqueContainer<std::string> &Methods);

	virtual void estimateInitialPriorParameters(const TSpecies &Data);

	void updateGamma(const TSpecies &Data);
	double calculateLL(const TSpecies &Data) const;

	// simulation
	void simulate(TSpecies &data, size_t N0, const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
	              const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods);
	void simulate(TSpecies &data);

	// get rho for tests
	size_t rho(size_t row, size_t col) const;
};

//--------------------------------
// TPoissonEpochModel
//--------------------------------

class TPoissonEpochModel : public TEpochModel {
protected:
	double _calculateLL(const std::vector<double> &sigma, const TDataSeries &lm, const std::vector<double> &p) const override;

public:
	TPoissonEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	                   stattools::TParameterTyped<TypeTime, 1> *TimesOfChange);
};

//--------------------------------
// TNegBinEpochModel
//--------------------------------

class TNegBinEpochModel : public TEpochModel {
protected:
	// parameters: gamma (in base class) and mu
	stattools::TParameterTyped<TypeMu, 2> *_mu = nullptr;

	double _calculateLL(const std::vector<double> &sigma, const TDataSeries &lm, const std::vector<double> &p) const override;

public:
	TNegBinEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	                  stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
	                  stattools::TParameterTyped<TypeMu, 2> *Mu);

	void initializeInferred(const TUniqueContainer<size_t> &Timepoints, const TUniqueContainer<std::string> &Locations,
	                        const TUniqueContainer<std::string> &Methods) override;
	void estimateInitialPriorParameters(const TSpecies &Data) override;
};

#endif // BIRP_TEPOCHMODEL_H
