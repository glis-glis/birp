//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "stattools/ParametersObservations/TObservationTyped.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include "TSpecies.h"
#include "coretools/Files/TFile.h"
#include "coretools/Files/TOutputFile.h"

class TBirpPrior : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2> {
private:
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::Storage;
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::UpdatedStorage;


protected:
	// data
	std::vector<TSpecies> _species;

	// dimensions
	TUniqueContainer<std::string> _locations;
	TUniqueContainer<std::string> _methods;
	TUniqueContainer<size_t> _timepoints;


	// simulate values
	void _simulateUnderPrior(Storage *Data) override;

	// check header and read file
	std::vector<std::string> _getHeader() const;
	void _checkHeader(const coretools::TInputFile &Infile);
	bool _checkHeaderPlacement(const coretools::TInputFile &Infile, size_t index, const std::string &Title);
	void _readFile(const std::string &Filename, bool hasCounts);
	coretools::TMultiDimensionalStorage<double, 2> _readDetectionProbFile();
	std::vector<size_t> _readN0();
	void _writeFile(const std::string &Filename);
	void _sortTimepoints();
	// write output file


public:
	TBirpPrior() = default;
	TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma);
	~TBirpPrior() override = default;

	[[nodiscard]] std::string name() const override;

	// temp
	void read();



	// initialization
	void initializeInferred() override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

    void simulate();
    // user: gives gamma, timepoints, times of change
    // counts and efforts: simulate


    // getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string>& locations() const;
    const TUniqueContainer<std::string>& methods() const;
    const TUniqueContainer<size_t>& timepoints() const;
    const std::vector<TSpecies>& species() const;

};

#endif // TBIRPPRIOR_H
