//
// Created by Liam on 23/08/22.
//

#include "TEpochModel.h"
#include "stattools/ParametersObservations/TIndexSampler.h"

using namespace coretools::instances;

//--------------------------------
// TEpochModel
//--------------------------------

TEpochModel::TEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                         stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
        : _gamma(Gamma), _timesOfChange(TimesOfChange) {

    // find out how many epochs there are
    if (parameters().parameterExists("timesOfChange")) {
        std::vector<std::string> tmp;
        parameters().fillParameterIntoContainer("timesOfChange", tmp, ',');
        _numEpochs = tmp.size() + 1;
    } else {
        _numEpochs = parameters().getParameterWithDefault("numEpochs", 1);
    }
    // default: do not infer times of change
    if (!parameters().parameterExists("inferTimesOfChange")) { _timesOfChange->getDefinition().update(false); }
}

void TEpochModel::initializeInferred(const TUniqueContainer<size_t> &Timepoints,
                                     const TUniqueContainer<std::string> & /*Locations*/,
                                     const std::vector<TMethods> &Methods) {
    _timesOfChange->initStorage({_numEpochs - 1}, {std::make_shared<coretools::TNamesIndices>()});
    _gamma->initStorage({_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});

    _fillRho(Timepoints);

    // check how many combinations of methods/locations we have - used to initialize vector of LL
    size_t numMethLoc = 0;
    for (auto m: Methods) {
        if (m.speciesExists(_species_id)){
            for (auto loc: m){
                numMethLoc = numMethLoc + 1;
            }
        }
    }
    _numMethLoc = numMethLoc;
}

void TEpochModel::estimateInitialPriorParameters(const std::vector<TMethods> &Methods) {
    // TODO: something more sophisticated?
    for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->set(g, 0.0); }

    // set _oldLL and _oldJP
    std::tie(_oldLL, _oldJP) = _calculateLLandJP(Methods);
}

std::vector<size_t> TEpochModel::_getEpochStarts(const TUniqueContainer<size_t> &Timepoints) const {
    // return start of each epoch, including for a hypothetical epoch after the last epoch
    assert(!Timepoints.empty());

    std::vector<size_t> epochStarts = {0};
    if (_timesOfChange->totalSize() != 0) {
        size_t maxTimeOfChange = _timesOfChange->max_value();
        size_t maxTime = std::max(Timepoints.back(), maxTimeOfChange + 1);
        for (size_t i = 0; i < _timesOfChange->totalSize(); i++) { epochStarts.push_back(_timesOfChange->value(i)); }
        epochStarts.push_back(maxTime);
    } else {
        epochStarts.push_back(*max_element(Timepoints.cbegin(), Timepoints.cend()));
    }
    return epochStarts;
}

void TEpochModel::_fillRho(const TUniqueContainer<size_t> &Timepoints) {
    std::vector<size_t> epochStarts = _getEpochStarts(Timepoints);
    _rho.zeros((int) Timepoints.size(), (int) _numEpochs);

    for (size_t t = 0; t < Timepoints.size(); ++t) {
        for (size_t e = 1; e < epochStarts.size(); ++e) {
            if (Timepoints[t] <= epochStarts[e - 1]) {
                _rho(t, e - 1) = 0;
            } else if (epochStarts[e - 1] < Timepoints[t] && Timepoints[t] < epochStarts[e]) {
                _rho(t, e - 1) = Timepoints[t] - epochStarts[e - 1];
            } else {
                _rho(t, e - 1) = epochStarts[e] - epochStarts[e - 1];
            }
        }
    }
}

double TEpochModel::_calculateLogPhi_k(size_t k) const {
    double sum = 0.0;
    for (size_t m = 0; m < _numEpochs; m++) { sum += (double) _rho(k, m) * _gamma->value(m); }
    return sum;
}

std::vector<double> TEpochModel::_calculatePhi() const {
    std::vector<double> phi(_rho.n_rows);
    for (size_t k = 0; k < _rho.n_rows; k++) { phi[k] = _calculateLogPhi_k(k); }
    coretools::standardizeZeroMean(phi);
    for (auto &p: phi) { p = exp(p); } // make logPhi to phi
    return phi;
}

arma::mat TEpochModel::_calcFisherInfo(const std::vector<double> &p, const std::vector<double> &sigma,
                                       TypeNu nu) const {
    arma::mat R(p.size(), _numEpochs);
    arma::mat RtP(_numEpochs, p.size());
    for (size_t k = 0; k < p.size(); k++) {
        for (size_t m = 0; m < _numEpochs; m++) {
            double denom = 0.0;
            double rhoSumSigma = 0.0;
            for (size_t l = 0; l < p.size(); l++) {
                denom += sigma[l];
                rhoSumSigma += ((double) _rho(k, m) - (double) _rho(l, m)) * sigma[l];
            }
            denom = denom * denom;
            R(k, m) = sigma[k] * rhoSumSigma / denom;
            RtP(m, k) = R(k, m) * 1.0 / p[k];
        }
    }
    // real actual matrix multiplication
    auto tmp = RtP * R;
    return nu * tmp;
}

std::vector<double> TEpochModel::calcSigma(const TLocations &loc, const std::vector<double> &phi) const {
    std::vector<double> sigma(loc.size());
    for (size_t k = 0; k < loc.size(); k++) {
        // calculate sigma and fill into p (will normalize later)
        sigma[k] = phi[k] * loc[k].effort();
    }
    return sigma;
}

std::pair<std::vector<double>, double> TEpochModel::_calculateLLandJP(const std::vector<TMethods> &Methods) const {
    std::vector<double> phi = _calculatePhi();
    std::vector<double> LL(_numMethLoc, 0.0);  // LL is per method/location combination
    arma::mat fisherInfo(_numEpochs, _numEpochs, arma::fill::zeros);
    size_t counter = 0;
    for (size_t m = 0; m < Methods.size(); m++) {
        if (Methods[m].speciesExists(_species_id)) {
            size_t index_in_counts = Methods[m].species_id_in_uniqueContainer(_species_id);
            for (size_t l = 0; l < Methods[m].size(); l++, counter++) {
                std::vector sigma = calcSigma(Methods[m][l], phi);
                std::vector<double> p = sigma;
                // now normalize by row sum!
                coretools::normalize(p);
                // add to Likelihood
                LL[counter] = _calculateLL(sigma, Methods[m][l], p, Methods[m].method_id(), index_in_counts);
                // add to JP
                fisherInfo += _calcFisherInfo(p, sigma, Methods[m][l].nu(_species_id));
            }
        }
    }
    double JP = log(arma::det(fisherInfo)); // sqrt drops out in hastings ratio
    return {LL, JP};
}

double TEpochModel::calculateLL(const std::vector<TMethods> &Methods) const {
    const auto [newLL, newJP] = _calculateLLandJP(Methods);
    return coretools::containerSum(newLL);
}

void TEpochModel::_updateGamma(const std::vector<TMethods> &Methods) {
    for (size_t m = 0; m < _numEpochs; m++) { // loop over epochs and update gamma epoch by epoch
        if (_gamma->updateSpecificIndex(m)) {
            const auto [newLL, newJP] = _calculateLLandJP(Methods);
            const double LLRatio = coretools::containerSum(newLL) - coretools::containerSum(_oldLL); // in log
            const double logPriorRatio = newJP - _oldJP;                                                   // in log

            const double logH = LLRatio + logPriorRatio;
            if (_gamma->acceptOrReject(logH, m)) {
                // accepted -> update LL and jeffrey's Prior
                _oldLL = newLL;
                _oldJP = newJP;
            }
        }
    }
}

void TEpochModel::update(const std::vector<TMethods>& Methods) { _updateGamma(Methods); }

void TEpochModel::_setSimulatedGammaToZero() {
    if (!_gamma->hasFixedInitialValue()) {
        for (size_t g = 0; g < _gamma->totalSize(); ++g) { _gamma->set(g, 0.0); }
    }
}

// getters for test functions
size_t TEpochModel::rho(size_t row, size_t col) const { return _rho(row, col); }

//--------------------------------
// TPoissonEpochModel
//--------------------------------

TPoissonEpochModel::TPoissonEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                                       stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
        : TEpochModel(Gamma, TimesOfChange) {};

double TPoissonEpochModel::_calculateLL(const std::vector<double> & /*sigma*/, const TLocations &loc,
                                        const std::vector<double> &p, size_t /*method_id*/, size_t index_in_counts) const {
    double LL = 0.0;
    for (size_t k = 0; k < loc.size(); k++) { LL += (double) loc[k].counts_per_species(index_in_counts) * log(p[k]); }
    return LL;
}

void TPoissonEpochModel::simulate(TMethods &data, size_t N0,
                                  const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
                                  const TUniqueContainer<std::string> &locations) {
    _setSimulatedGammaToZero();
    data.simulatePoisson(N0, _calculatePhi(), detectionProbs, locations);
}

void TPoissonEpochModel::simulate(TMethods &data) {
    _setSimulatedGammaToZero();
    data.simulatePoisson(_calculatePhi());
}

//--------------------------------
// TNegBinEpochModel
//--------------------------------

TNegBinEpochModel::TNegBinEpochModel(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                                     stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                                     stattools::TParameterTyped<TypeMu, 2> *Mu,
                                     stattools::TParameterTyped<TypeLogA, 1> *LogA)
        : TEpochModel(Gamma, TimesOfChange) {
    _mu = Mu;
    _logA = LogA;

    // set max of logA to something reasonable: otherwise we have numeric issues with gammaLog function
    TypeLogA::setMax(10.0);
}

void TNegBinEpochModel::initializeInferred(const TUniqueContainer<size_t> &Timepoints,
                                           const TUniqueContainer<std::string> &Locations,
                                           const std::vector<TMethods> &Methods) {
    // initialize gamma and timepoints
    TEpochModel::initializeInferred(Timepoints, Locations, Methods);

    // initialize mu
    auto locationNames = std::make_shared<coretools::TNamesStrings>();
    auto methodNames = std::make_shared<coretools::TNamesStrings>();
    for (size_t j = 0; j < Locations.size(); ++j) { locationNames->addName({Locations[j]}); }
    for (size_t i = 0; i < Methods.size(); ++i) { methodNames->addName({Methods[i].name()}); }
    _mu->initStorage({Methods.size(), Locations.size()}, {methodNames, locationNames});

    // initialize a: single element
    _logA->initStorage({1});
}

double TNegBinEpochModel::_calculateLL(const std::vector<double> &sigma, const TLocations &loc,
                                       const std::vector<double> & /*p*/, size_t method_id,
                                       size_t index_in_counts) const {
    double sum_a = 0.0;
    double term2 = 0.0;
    TypeNu n = loc.nu(index_in_counts);
    const double a = exp(_logA->value());
    const size_t ij = _mu->getIndex({method_id, loc.location_id()});
    for (size_t k = 0; k < loc.size(); k++) {
        double alpha_k = a * _mu->value(ij) * sigma[k];
        sum_a += alpha_k; // summing up alphas
        term2 += coretools::gammaLog((double) loc[k].counts_per_species(index_in_counts) + alpha_k) -
                 coretools::gammaLog(alpha_k); // x_k + 1 drops out in hastings ratio
    }
    return coretools::gammaLog(sum_a) - coretools::gammaLog((double) n + sum_a) + term2;  // (n+1) drops out
}


std::pair<size_t, size_t> TNegBinEpochModel::_getMethAndLocIndex(size_t mu_index, const std::vector<TMethods>& Methods){
    size_t meth_id = 0;
    size_t loc_id = 0;
    size_t counter = 0;
    for (const auto& method: Methods){
        if(method.speciesExists(_species_id)) {
            if (counter + method.size() > mu_index) {
                meth_id = method.method_id();
                loc_id = mu_index - counter;
                return {meth_id, loc_id};
            }
        }
        counter += method.size();
    }
    DEVERROR("Error in TNegBinEpochModel::_getMethAndLocIndex(): did not find method and location index for species - this should not happen.");
}

void TNegBinEpochModel::_updateMu(const std::vector<TMethods> &Methods) {
    stattools::TIndexSampler sampler(_numMethLoc);
    sampler.sampleIndices();

    for (size_t c = 0; c < sampler.length(); ++c) {
        // pick two indices i to update
        auto [c1, c2] = sampler.getIndexPair(c);

        auto [m_id_c1, l_id_c1] = _getMethAndLocIndex(c1, Methods);
        auto [m_id_c2, l_id_c2] = _getMethAndLocIndex(c2, Methods);
        const auto &loc1 = Methods[m_id_c1][l_id_c1];
        const auto &loc2 = Methods[m_id_c2][l_id_c2];

        const size_t ij1 = _mu->getIndex({Methods[m_id_c1].method_id(), loc1.location_id()});
        const size_t ij2 = _mu->getIndex({Methods[m_id_c2].method_id(), loc2.location_id()});

        _mu->updatePairWiseSumOneVector(ij1, ij2);

        std::vector<double> phi = _calculatePhi();
        std::vector<double> sigma1 = calcSigma(loc1, phi);
        std::vector<double> sigma2 = calcSigma(loc2, phi);

        size_t index_in_counts1 = Methods[m_id_c1].species_id_in_uniqueContainer(_species_id);
        const auto newLL1 = _calculateLL(sigma1, loc1, {}, Methods[m_id_c1].method_id(), index_in_counts1);
        size_t index_in_counts2 = Methods[m_id_c1].species_id_in_uniqueContainer(_species_id);
        const auto newLL2 = _calculateLL(sigma2, loc2, {}, Methods[m_id_c2].method_id(), index_in_counts2);

        const double LLRatio = newLL1 + newLL2 - _oldLL[c1] - _oldLL[c2]; // in log
        const double logH = LLRatio + _mu->getLogPriorRatio(ij1) + _mu->getLogPriorRatio(ij2);

        coretools::TRange range(ij1, ij2);
        if (_mu->acceptOrReject(logH, range, ij1)) {
            // accepted -> update LL
            _oldLL[c1] = newLL1;
            _oldLL[c2] = newLL2;
        }
    }
}

void TNegBinEpochModel::_updateLogA(const std::vector<TMethods> &Methods) {
    if (_logA->update()) {
        std::vector<double> phi = _calculatePhi();
        std::vector<double> newLL(_numMethLoc);
        for (size_t ij = 0; ij < _numMethLoc; ++ij) {
            auto [m_id, l_id] = _getMethAndLocIndex(ij, Methods);
            const auto &loc = Methods[m_id][l_id];
            std::vector sigma = calcSigma(loc, phi);
            size_t index_in_counts = Methods[m_id].species_id_in_uniqueContainer(_species_id);
            newLL[ij] = _calculateLL(sigma, loc, {}, Methods[m_id].method_id(), index_in_counts);
        }
        const double LLRatio = coretools::containerSum(newLL) - coretools::containerSum(_oldLL);
        const double logH = LLRatio + _logA->getLogPriorRatio();
        if (_logA->acceptOrReject(logH, 0)) {
            // accepted -> update LL
            _oldLL = newLL;
        }
    }
}

void TNegBinEpochModel::update(const std::vector<TMethods> &Methods) {
    if (_gamma->isUpdated()) { _updateGamma(Methods); }
    if (_mu->isUpdated()) { _updateMu(Methods); }
    if (_logA->isUpdated()) { _updateLogA(Methods); }
}

void TNegBinEpochModel::_estimateInitialMu(const std::vector<TMethods> &Methods) {
    // initialize mu: relative nu's
    std::vector<double> nus(_numMethLoc);
    for (size_t c = 0; c < _numMethLoc; ++c) {
        auto [m_id, l_id] = _getMethAndLocIndex(c, Methods);
        nus[c] = std::max(coretools::UnsignedInt<uint32_t>(1),
                          Methods[m_id][l_id].nu(_species_id)); // do not allow 0: problem for Dirichlet
    }
    coretools::normalize(nus);

    // initialize to very small value: not all i-j-combinations exists
    for (size_t ij = 0; ij < _mu->totalSize(); ++ij) { _mu->setInitial(ij, 1e-10); }
    for (size_t c = 0; c <  _numMethLoc; ++c) {
        auto [m_id, l_id] = _getMethAndLocIndex(c, Methods);
        const size_t ij = _mu->getIndex({Methods[m_id].method_id(), Methods[m_id][l_id].location_id()});
        _mu->setInitial(ij, nus[c]);
    }
    _mu->normalizeToSum1AfterUpdate();
}



void TNegBinEpochModel::estimateInitialPriorParameters(const std::vector<TMethods> &Methods) {
    // initialize mu
    _estimateInitialMu(Methods);
    _logA->setInitial(0.0); // log(1) = 0

    TEpochModel::estimateInitialPriorParameters(Methods);
}


double TNegBinEpochModel::_simulate(TMethods &Methods) {
    _setSimulatedGammaToZero();
    if (!_logA->hasFixedInitialValue()) { _logA->set(log(Methods.size())); }
    // calculate a = exp(log(a))
    const double a = exp(_logA->value());
    return a;
}


void TNegBinEpochModel::simulate(TMethods &data, size_t N0,
                                 const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
                                 const TUniqueContainer<std::string> &locations) {
    double A = _simulate(data);
    std::vector<double> phi = _calculatePhi();
    auto mus = data.simulateNB(N0, phi, detectionProbs, locations, A); // simulate mus based on detectionProb and a
    coretools::normalize(mus);
    for (size_t c = 0; c < data.size(); ++c) {
        const size_t ij = _mu->getIndex({data.method_id(), data[c].location_id()});
        _mu->set(ij, mus[c]); // store those mus here
    }
}


void TNegBinEpochModel::simulate(TMethods &data) {
    double a = _simulate(data);
    // collect mu from parameter (only relevant method-location combinations)
    _mu->normalizeToSum1AfterUpdate();
    std::vector<double> mu(data.size());
    for (size_t c = 0; c < data.size(); ++c) {
        const size_t ij = _mu->getIndex({data.method_id(), data[c].location_id()});
        mu[c] = _mu->value(ij);
    }
    data.simulateNB(_calculatePhi(), a, mu, true);

}

// TODO
// automatically remove data sets with less than two time points with effort > 0?
