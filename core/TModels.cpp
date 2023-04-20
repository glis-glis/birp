//
// Created by Liam on 23/08/22.
//
#include "TModels.h"
#include "TDistributions.h"

//---------------------------------------
// TModelBase
//---------------------------------------

TModelBase::TModelBase(size_t SpeciesIDinUniqueContainer, const TData &Data,
                       const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn)
    : _speciesIDinUniqueContainer(SpeciesIDinUniqueContainer) {

	// create distribution
	if (NegBinomAddOn) { // negative binomial
		_distr = std::make_unique<TNegBinDistribution>(_speciesIDinUniqueContainer, &NegBinomAddOn->mu,
		                                               &NegBinomAddOn->logA);
	} else { // Poisson
		_distr = std::make_unique<TPoissonDistribution>(_speciesIDinUniqueContainer);
	}

	// LL is per method/location combination
	_curLL.resize(Data.numMethLoc(_speciesIDinUniqueContainer), 0.0);
	_tryLL.resize(Data.numMethLoc(_speciesIDinUniqueContainer), 0.0);
}

void TModelBase::_initializeDistribution(const TData &Data, const TUniqueContainer<std::string> &Locations) {
	_distr->initialize(Data, Locations);
}

void TModelBase::_estimateInitialParameters(const TData &Data) { _distr->estimateInitialParameters(Data); }

void TModelBase::_updateDistributionParameters(const TData &Data) { _distr->update(Data, *this); };

void TModelBase::_simulateDistributionParameters(
    TData &Data, size_t N0,
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>> &detectionProbs,
    const TUniqueContainer<std::string> &locations, bool SimulateCovariates) {

	_distr->simulate(Data, *this, N0, detectionProbs, locations, SimulateCovariates);
}

std::vector<double> TModelBase::calculateSigma(const TLocations &loc) const {
	std::vector<double> sigma(loc.size());
	for (size_t k = 0; k < loc.size(); k++) {
		// calculate sigma and fill into p (will normalize later)
		sigma[k] = phi(loc.location_id(), loc[k].timepoint_id()) * loc[k].effort();
	}
	return sigma;
}

arma::mat TModelBase::_calcFisherInfo(const TGamma &Gamma, const std::vector<double> &p,
                                      const std::vector<double> &sigma, TypeNu nu) {
	arma::mat R(p.size(), Gamma.numEpochs());
	arma::mat RtP(Gamma.numEpochs(), p.size());
	for (size_t k = 0; k < p.size(); k++) {
		for (size_t m = 0; m < Gamma.numEpochs(); m++) {
			double denom       = 0.0;
			double rhoSumSigma = 0.0;
			for (size_t l = 0; l < p.size(); l++) {
				denom += sigma[l];
				rhoSumSigma += ((double)Gamma.rho(k, m) - (double)Gamma.rho(l, m)) * sigma[l];
			}
			denom     = denom * denom;
			R(k, m)   = sigma[k] * rhoSumSigma / denom;
			RtP(m, k) = R(k, m) * 1.0 / p[k];
		}
	}
	// real actual matrix multiplication
	auto tmp = RtP * R;
	return nu * tmp;
}

double TModelBase::calculateJeffreyPrior(const TData &Data, const TGamma &Gamma) const {
	// used in update of gamma of stochastic model
	arma::mat fisherInfo(Gamma.numEpochs(), Gamma.numEpochs(), arma::fill::zeros);
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = calculateSigma(method[l]);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to JP
				fisherInfo += _calcFisherInfo(Gamma, p, sigma, method[l].nu(index_in_counts));
			}
		}
	}
	return log(arma::det(fisherInfo)); // sqrt drops out in hastings ratio
}

void TModelBase::_initializeLL(const TData &Data) {
	// loop over all datapoints (methods - locations - timepoints) of a single species
	// used to initialize curLL
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = calculateSigma(method[l]);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to Likelihood
				_curLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
			}
		}
	}
};

double TModelBase::calculateLLRatio(const TData &Data) {
	// loop over all datapoints (methods - locations - timepoints) of a single species
	// used in update of gamma of deterministic model and sigma in stochastic model
	size_t ix_linNumMethLoc = 0;
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
			for (size_t l = 0; l < method.size(); l++, ix_linNumMethLoc++) { // loop over locations
				std::vector sigma     = calculateSigma(method[l]);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to Likelihood
				_tryLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc);
			}
		}
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
};

double TModelBase::calculateLLRatio_perLocation(size_t j, const TData &Data) {
	// loop over all methods and timepoints of a particular location of a single species
	// used in update of phi in stochastic model
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			size_t index_in_counts  = method.getIndexInCounts(_speciesIDinUniqueContainer);
			std::vector sigma       = calculateSigma(method[j]);
			std::vector<double> p   = sigma;
			// now normalize by row sum!
			coretools::normalize(p);
			// add to Likelihood
			_tryLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
		}
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
};

double TModelBase::calculateLLRatio_perMethod(size_t i, const TData &Data) {
	// loop over all locations and timepoints of a particular method of a single species
	const auto &method = Data[i];

	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
		return coretools::containerSum(_curLL);
	}                                            // chame gd hei ga

	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		const auto &loc        = method[j];
		size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
		std::vector sigma      = calculateSigma(loc);
		std::vector<double> p  = sigma;
		// now normalize by row sum!
		coretools::normalize(p);
		// calculate Likelihood
		size_t ix_linNumMethLoc  = Data.i_j_to_linear(_speciesIDinUniqueContainer, method.method_id(),
		                                              j); // method_id should always be the same as i!
		_tryLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc);
	}
	return coretools::containerSum(_tryLL);
}

void TModelBase::swapTryCur() {
	// used when all datapoints are affected (update of deterministic gamma or sigma)
	std::swap(_tryLL, _curLL);
}

void TModelBase::swapTryCur_perLocation(size_t j, const TData &Data) {
	// used when a location has changed
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &method = Data[i];
		if (method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			std::swap(_tryLL[ix_linNumMethLoc], _curLL[ix_linNumMethLoc]);
		}
	}
}

void TModelBase::swapTryCur_perMethod(size_t i, const TData &Data) {
	// used when a method has changed
	const auto &method = Data[i];
	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }
	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, method.method_id(), j);
		std::swap(_tryLL[ix_linNumMethLoc], _curLL[ix_linNumMethLoc]);
	}
}

double TModelBase::curLL(size_t ij) const { return _curLL[ij]; }
double TModelBase::sumCurLL() const { return coretools::containerSum(_curLL); }
void TModelBase::setCurLL(size_t ij, double Value) { _curLL[ij] = Value; }
void TModelBase::setCurLL(const std::vector<double> &Values) { _curLL = Values; }

// TODO: issue warning if data sets with less than two time points with effort > 0

//------------------------------------------
// TDeterministic
//------------------------------------------

TDeterministic::TDeterministic(size_t SpeciesIDinUniqueContainer, const TData &Data,
                               const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
                               stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                               stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
    : TModelBase(SpeciesIDinUniqueContainer, Data, NegBinomAddOn), _gamma(Gamma, TimesOfChange) {}

void TDeterministic::initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
                                const TUniqueContainer<std::size_t> &TimePoints) {
	_initializeDistribution(Data, Locations);
	_gamma.initialize(TimePoints);

	_newPhi.resize(_gamma.numTimepoints());
	_oldPhi = _newPhi;
}

void TDeterministic::estimateInitial(const TData &Data) {
	_estimateInitialParameters(Data);
	_gamma.estimateInitialGamma();
	_fillPhi();
	_oldPhi = _newPhi;
}

void TDeterministic::initializeTempVariables(const TData &Data) {
	_gamma.initializeJP(Data, *this);
	_initializeLL(Data);
}

void TDeterministic::_fillPhi() { _newPhi = _gamma.phi(); }

TypePhi TDeterministic::phi(size_t /*location_id*/, size_t timepoint_id) const { return _newPhi[timepoint_id]; }

void TDeterministic::update(const TData &Data) {
	// update distribution parameters
	_updateDistributionParameters(Data);

	for (size_t m = 0; m < _gamma.numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma.updateSpecificIndex(m)) {
			_fillPhi();
			const auto LLRatio         = calculateLLRatio(Data);
			const double logPriorRatio = _gamma.getLogPriorRatio(Data, *this); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma.acceptOrReject(logH, m)) {
				// accepted
				swapTryCur();
			} else {
				// reject
				_newPhi = _oldPhi;
			}
		}
	}
}

void TDeterministic::simulate(
    TData &Data, size_t N0,
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>> &detectionProbs,
    const TUniqueContainer<std::string> &locations, bool SimulateCovariates) {

	_gamma.setSimulatedGammaToZero();
	_fillPhi();
	_simulateDistributionParameters(Data, N0, detectionProbs, locations, SimulateCovariates);
};

//------------------------------------------
// TStochastic
//------------------------------------------

TStochastic::TStochastic(size_t SpeciesIDinUniqueContainer, const TData &Data,
                         const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
                         stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                         stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                         stattools::TParameterTyped<TypePhi, 2> *Phi, stattools::TParameterTyped<TypeSigma, 1> *Sigma)
    : TModelBase(SpeciesIDinUniqueContainer, Data, NegBinomAddOn), _prior(Phi, Sigma, Gamma, TimesOfChange), _phi(Phi),
      _sigma(Sigma) {}

void TStochastic::initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
                             const TUniqueContainer<std::size_t> &TimePoints) {
	_initializeDistribution(Data, Locations);
	_prior.initialize(TimePoints);

	// create name classes for locations
	auto locNames = std::make_shared<coretools::TNamesStrings>(Locations.vec());
	std::vector<std::string> tmp(TimePoints.size());
	for (size_t t = 0; t < TimePoints.size(); ++t) { tmp[t] = coretools::str::toString(TimePoints[t]); }
	auto timeNames = std::make_shared<coretools::TNamesStrings>(tmp);

	_phi->initStorage({Locations.size(), TimePoints.size()}, {locNames, timeNames});
	_sigma->initStorage();
}

void TStochastic::_setPhiFromDeterministic() {
	auto deterministicPhi = _prior.phi();
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 0; k < _phi->dimensions()[1]; ++k) { // loop over timepoints
			const auto index = _phi->getIndex({j, k});       // get linear index
			_phi->setInitial(index, deterministicPhi[k]);
		}
	}
}

void TStochastic::estimateInitial(const TData &Data) {
	_estimateInitialParameters(Data);
	_prior.estimateInitialGamma();

	// initialize phi based on deterministic phi
	_setPhiFromDeterministic();

	// set smallish sigma
	_sigma->setInitial(0.1);
}

void TStochastic::initializeTempVariables(const TData &Data) {
	_prior.initializeTempVariables(Data, *this);
	_initializeLL(Data);
}

void TStochastic::_updatePhi(size_t j, size_t k, const TData &Data) {
	// propose a new value
	const size_t ix = _phi->getIndex({j, k});
	_phi->updateSpecificIndex(ix);

	const auto logPriorRatio = _prior.calculateLogPriorRatio(j, k) + _prior.calculateLogPriorRatio(j, k + 1);
	// only need to calculate LL for this location j and timepoint k, but across all methods (if they have data for j
	// and k)
	const auto LLRatio       = calculateLLRatio_perLocation(j, Data);
	const double logH        = LLRatio + logPriorRatio;
	if (_phi->acceptOrReject(logH, ix)) {
		_prior.swapTryCur(j, k);
		_prior.swapTryCur(j, k + 1);
		swapTryCur_perLocation(j, Data);
	}
}

void TStochastic::_updatePhi(const TData &Data) {
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _phi[j].dimensions()[1]; ++k) {
			// loop over all timepoints (ignore first, since this is always 1)
			_updatePhi(j, k, Data);
		}
	}
}

void TStochastic::_updateSigma() {
	// propose a new value
	if (!_sigma->update()) { return; }

	double LLRatio = 0.0;
	for (size_t j = 0; j < _phi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _phi[j].dimensions()[1]; ++k) { LLRatio += _prior.calculateLogPriorRatio(j, k); }
	}

	const double logH = LLRatio + _sigma->getLogPriorRatio();
	if (_sigma->acceptOrReject(logH, 0)) {
		// accepted: update transition probabilities
		_prior.swapTryCur();
		swapTryCur();
	}
}

void TStochastic::update(const TData &Data) {
	// update distribution parameters
	_updateDistributionParameters(Data);

	// update gamma
	_prior.updateGamma(Data, *this);

	if (_phi->isUpdated()) { _updatePhi(Data); }
	if (_sigma->isUpdated()) { _updateSigma(); }
}

void TStochastic::simulate(
    TData &Data, size_t N0,
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>> &detectionProbs,
    const TUniqueContainer<std::string> &locations, bool SimulateCovariates) {

	// simulate gamma
	_prior.setSimulatedGammaToZero();

	_setPhiFromDeterministic();

	_simulateDistributionParameters(Data, N0, detectionProbs, locations, SimulateCovariates);
};

coretools::StrictlyPositive<double> TStochastic::phi(size_t j, size_t k) const {
	return _phi->value(_phi->getIndex({j, k}));
}
