
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

    // temporary variables
	std::vector<double> _phi;

    size_t _numEpochs = 0;
    coretools::TMatrix<size_t> _rho;
    std::vector<TypeCounts> _nu; // sums of counts at one location and method (ij) across all time points (k)

    // other stuff
    std::string _filename;

    // likelihood
    [[nodiscard]] double _calculateLogPhi_k(size_t k) const;
    [[nodiscard]] std::vector<double> _calculatePhi(size_t numTimepoints) const;
    double _calculateSigma_ik(size_t i, size_t k, double Phi_k) const;
    coretools::TMatrix<double> _calculateP(size_t numTimepoints, size_t numMethods, const TDataSeries &Data) const;

	void _fillRho(const TUniqueContainer<size_t> &Timepoints);


public:
    TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1>* TimesOfChange, const std::string &Filename);
    ~TEpochModel() = default;

	void initialize();

    // initialization
    [[nodiscard]] std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const; //return start of each epoch, including for a hypothetical epoch after the last epoch

    // static void _readTimesOfChange(); already done in readFile?

    double calculateLL(const TDataSeries &Data);

	// getter for rho for tests
	size_t rho(size_t row, size_t col) const;
    // const  std::vector<size_t>& epochStarts() const;


};

#endif // BIRP_TEPOCHMODEL_H

