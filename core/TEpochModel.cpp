//
// Created by Liam on 23/08/22.
//

#include "TEpochModel.h"
#include "TData.h"

using namespace coretools::instances;

TEpochModel::TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                         stattools::TObservationTyped<TypeEffort, 2> *Effort, const std::string &Filename)
        : _gamma(Gamma), _effort(Effort), _filename(Filename) {
    addPriorParameter({_gamma, _effort});
}

std::vector<size_t> TEpochModel::_getEpochStarts(const std::vector<size_t> &TimesOfChange,
                                                const std::vector<size_t> &Times) {
    assert(!Times.empty());
    std::vector<size_t> epochStarts = {0};
    if (!TimesOfChange.empty()) {
        size_t maxTimeOfChange = *max_element(TimesOfChange.begin(), TimesOfChange.end());
        size_t maxTime         = std::max(Times.back(), maxTimeOfChange + 1);
        //
        epochStarts.insert(epochStarts.end(), TimesOfChange.begin(), TimesOfChange.end());
        epochStarts.push_back(maxTime);
    } else {
        epochStarts.push_back(*max_element(Times.begin(), Times.end()));
    }
    return epochStarts;
}


void TEpochModel::_fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
    std::vector<size_t> epochStarts = _getEpochStarts(TimesOfChange, Times);
    _numEpochs                      = TimesOfChange.size() + 1;
    _rho.zeros((int)Times.size(), (int)_numEpochs);
    // fill rhos
    for (size_t t = 0; t < Times.size(); ++t) {
        for (size_t e = 1; e < epochStarts.size(); ++e) {
            if (Times[t] <= epochStarts[e - 1]) {
                _rho(t, e - 1) = 0;
            } else if (epochStarts[e - 1] < Times[t] && Times[t] < epochStarts[e]) {
                _rho(t, e - 1) = Times[t] - epochStarts[e - 1];
            } else {
                _rho(t, e - 1) = epochStarts[e] - epochStarts[e - 1];
            }
        }
    }
}

void TEpochModel::_readTimesOfChange() {
    // read vector from command line
    std::vector<size_t> timesOfChange;
    parameters().fillParameterIntoContainer("timesOfChange", timesOfChange, ',');
    logfile().list("Will consider the following times of change: ",
                   coretools::str::concatenateString(timesOfChange, ","));

    // convert time points (strings, because they were read from file) into size_t
    const size_t numTimePoints = _effort->getDimensionName(1)->size();
    std::vector<size_t> times(numTimePoints);
    for (size_t t = 0; t < numTimePoints; t++) {
        times[t] = coretools::str::convertString<size_t>((*_effort->getDimensionName(1))[t]);
    }

    // now fill rho!
    _fillRho(timesOfChange, times);
}


double TEpochModel::_calculateLogPhi_k(size_t k) const {
    double sum = 0.0;
    for (size_t m = 0; m < _numEpochs; m++) { sum += (double)_rho(k, m) * _gamma->value(m); }
    return sum;
}

std::vector<double> TEpochModel::_calculatePhi() const {
    std::vector<double> logPhi(_numTimePoints);
    for (size_t k = 0; k < _numTimePoints; k++) { logPhi[k] = _calculateLogPhi_k(k); }
    coretools::standardizeZeroMean(logPhi);
    for (auto &p : logPhi) { p = exp(p); }
    return logPhi; // is not logPhi anymore due to exp() on line above
}


double TEpochModel::_calculateSigma_ik(size_t i, size_t k, double Phi_k) const {
    return Phi_k * _effort->value({i, k});
}

coretools::TMatrix<double> TEpochModel::_calculateP() const {
    std::vector<double> phi = _calculatePhi();

    coretools::TMatrix<double> p((int)_numMethods, (int)_numTimePoints);
    for (size_t i = 0; i < _numMethods; i++) {
        double sum = 0.0;
        for (size_t k = 0; k < _numTimePoints; k++) {
            // calculate sigma and fill into p (will normalize later)
            const double sigma = _calculateSigma_ik(i, k, phi[k]);
            p(i, k)            = sigma;
            sum += sigma;
        }
        // now normalize by row sum!
        for (size_t k = 0; k < _numTimePoints; k++) { p(i, k) /= sum; }
    }
    return p;
}

