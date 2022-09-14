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

	// data
	std::vector<TSpecies> _species;

	// dimensions
	TUniqueContainer<std::string> _locations;
	TUniqueContainer<std::string> _methods;
	TUniqueContainer<size_t> _timepoints;

protected:
	// simulate values
	void _simulateUnderPrior(Storage *Data) override;

	// check header and read file
	void _checkHeader(const coretools::TInputFile &Infile);
	bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title);
	void _readFile(const std::string &Filename);

public:
	TBirpPrior() = default;
	TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma);
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

	// getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string>& locations() const;
    const TUniqueContainer<std::string>& methods() const;
    const TUniqueContainer<size_t>& timepoints() const;
    const std::vector<TSpecies>& species() const;

};

#endif // TBIRPPRIOR_H
