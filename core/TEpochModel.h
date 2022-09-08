/*
//
// Created by Liam on 23/08/22.
//

#ifndef BIRP_TEPOCHMODEL_H
#define BIRP_TEPOCHMODEL_H

#include "BirpTypes.h"
#include "TMatrix.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"
#include "TSpecies.h"
#include <vector>


class TEpochModel {
private:
	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1>* _gamma = nullptr;
	stattools::TParameterTyped<TypeTime, 1>* _timesOfChange  = nullptr;
	std::vector<TypeTime> _observationTimes;

    // temporary variables
	std::vector<double> _phi;

    size_t _numEpochs = 0;
    coretools::TMatrix<size_t> _rho;
    std::vector<TypeCounts> _nu; // sums of counts at one location and method (ij) across all time points (k)

    // other stuff
    std::string _filename;

    // likelihood
    double _calculateLogPhi_k(size_t k) const;
    std::vector<double> _calculatePhi() const;
    double _calculateSigma_ik(size_t i, size_t k, double Phi_k) const;
    coretools::TMatrix<double> _calculateP() const;

public:
    TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TObservationTyped<TypeEffort, 2> *Effort,
    const std::string &Filename);
    ~TEpochModel() = default;

    // initialization
    std::vector<size_t> _getEpochStarts(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times); //return start of each epoch, including for a hypothetical epoch after the last epoch
    void _fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times);
    static void _readTimesOfChange();

    double calculateLL(const TData & Data);

};

#endif // BIRP_TEPOCHMODEL_H
*/
