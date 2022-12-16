
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
#include <vector>

class TEpochModel {
protected:
	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1>* _gamma = nullptr;
	stattools::TParameterTyped<TypeTime, 1>* _timesOfChange  = nullptr;

	double _oldLL = 0.0;

    size_t _numEpochs = 0;
    coretools::TMatrix<size_t> _rho;
    std::vector<TypeCounts> _nu; // sums of counts at one location and method (ij) across all time points (k)

    // other stuff
    std::string _filename;

    // likelihood
    [[nodiscard]] double _calculateLogPhi_k(size_t k) const;
    [[nodiscard]] std::vector<double> _calculatePhi() const;
    double _calculateLL(const TSpecies &Data) const;
	void _fillRho(const TUniqueContainer<size_t> &Timepoints);


public:
    TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1>* TimesOfChange, const std::string &Filename);
    ~TEpochModel() = default;

	void initialize();

    // initialization
    [[nodiscard]] std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const; //return start of each epoch, including for a hypothetical epoch after the last epoch

    // static void _readTimesOfChange(); already done in readFile?

    double calculateLL(const TDataSeries &Data);

	// getter for rho and numEpochs for tests
	size_t rho(size_t row, size_t col) const;



};

#endif // BIRP_TEPOCHMODEL_H

