//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "TMatrix.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"
#include "TSpecies.h"




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
	coretools::TMatrix<size_t> _rho;
	std::vector<TypeCounts> _nu; // sums of counts at one location and method (ij) across all time points (k)

	// names
	TUniqueContainer<std::string> _locations;
	TUniqueContainer<std::string> _methods;
	TUniqueContainer<size_t> _timepoints;
	std::vector<TSpecies> _species;

	// temporary variables, for speed
	double _oldLL           = 0.0;
	double _oldLogJeffPrior = 0.0;

	// other stuff
	std::string _filename;

protected:
	// simulate values
	void _simulateUnderPrior(Storage *Data) override;

	// likelihood
	double _calculateLL(const coretools::TMatrix<double> & p) const;
	double _calculateLL() const;

	// Jeffrey's prior
	double _calculateLogJeffreysPrior(const coretools::TMatrix<double> & p) const;

	// MCMC updates
	void _evaluateUpdateGamma(size_t m);

	// check header and read file
	void _checkHeader(const coretools::TInputFile &Infile);
	bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title);
	void _readFile();

public:
	TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TObservationTyped<TypeEffort, 2> *Effort,
	           const std::string &Filename);
	~TBirpPrior() override = default;

	[[nodiscard]] std::string name() const override;

	// initialization
	void initializeInferred() override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

	// getters
	double rho(size_t Time, size_t Epoch) const;
	coretools::TMatrix<size_t> rhoMat() const;



};



#endif // TBIRPPRIOR_H
