//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "TMatrix.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"

class TBirpPrior : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2> {
private:
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::Storage;
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::UpdatedStorage;

	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1> *_gamma     = nullptr;
	stattools::TObservationTyped<TypeEffort, 2> *_effort = nullptr;

	// main dimensions
	size_t _numMethods    = 0;
	size_t _numTimePoints = 0;

	// epochs, rho
	size_t _numEpochs = 0;
	coretools::TMatrix<double> _rho;

	// temporary variables, for speed
	double _oldLL           = 0.0;
	double _oldLogJeffPrior = 0.0;

	// other stuff
	std::string _filename;

protected:
	// simulate values
	void _simulateUnderPrior(Storage *Data) override;

	// initialization
	void _fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times);
	void _readTimesOfChange();

	// likelihood
	double _calculatePhi_k(size_t k) const;
	std::vector<double> _calculatePhi() const;
	double _calculateSigma_ik(size_t i, size_t k, double Phi_k) const;
	coretools::TMatrix<double> _calculateP() const;
	double _calculateLL(Storage *Counts) const;

	// Jeffrey's prior
	double _calculateLogJeffreysPrior() const;

	// MCMC updates
	void _evaluateUpdateGamma(size_t m, Storage *Counts);

public:
	TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TObservationTyped<TypeEffort, 2> *Effort,
	           const std::string &Filename);
	~TBirpPrior() override = default;

	[[nodiscard]] std::string name() const override;

	// initialization
	void initializeStorageOfPriorParameters() override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

	// getters
	double rho(size_t Time, size_t Epoch);
};

#endif // TBIRPPRIOR_H
