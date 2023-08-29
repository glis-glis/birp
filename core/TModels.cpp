//
// Created by Liam on 23/08/22.
//
#include "TModels.h"
#include "TDistributions.h"
#include "coretools/devtools.h"

//---------------------------------------
// TModelBase
//---------------------------------------

TModelBase::TModelBase(size_t SpeciesIDinUniqueContainer, const TData &Data,
                       const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn)
    : _speciesIDinUniqueContainer(SpeciesIDinUniqueContainer) {

	// create distribution
	if (NegBinomAddOn) { // negative binomial
		_distr = std::make_unique<TNegBinDistribution>(_speciesIDinUniqueContainer, &NegBinomAddOn->muOrN,
		                                               &NegBinomAddOn->b);
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

void TModelBase::_simulateDistributionParameters(TData &Data, const std::vector<double> &n_bar,
                                                 bool SimulateCovariates) {
	_distr->simulate(Data, *this, n_bar, SimulateCovariates);
}

std::vector<double> TModelBase::calculateSigma(const TLocations &loc) const {
	std::vector<double> sigma(loc.size());
	for (size_t k = 0; k < loc.size(); k++) {
		// calculate sigma and fill into p (will normalize later)
		sigma[k] = phi(loc.location_id(), loc[k].timepoint_id()) * loc[k].effort();
	}
	return sigma;
}

std::vector<double> TModelBase::calculateSigmaFromDeterministicPhi(const TLocations &loc) const {
	std::vector<double> sigma(loc.size());
	for (size_t k = 0; k < loc.size(); k++) {
		// calculate sigma and fill into p (will normalize later)
		sigma[k] = deterministicPhi(loc[k].timepoint_id()) * loc[k].effort();
	}
	return sigma;
}

arma::mat TModelBase::_calcFisherInfo(const TGamma &Gamma, const std::vector<double> &p,
                                      const std::vector<double> &sigma, TypeNu nu, const TLocations &loc) {
	// size of p and sigma = number of timepoints of current location
	arma::mat R(p.size(), Gamma.numEpochs());
	arma::mat RtP(Gamma.numEpochs(), p.size());
	for (size_t k = 0; k < p.size(); k++) {
		for (size_t m = 0; m < Gamma.numEpochs(); m++) {
			double denom       = 0.0;
			double rhoSumSigma = 0.0;
			for (size_t l = 0; l < p.size(); l++) {
				denom += sigma[l];
				const auto diffRho =
				    (double)Gamma.rho(loc[k].timepoint_id(), m) - (double)Gamma.rho(loc[l].timepoint_id(), m);
				rhoSumSigma += diffRho * sigma[l];
			}
			denom     = denom * denom;
			R(k, m)   = sigma[k] * rhoSumSigma / denom;
			RtP(m, k) = R(k, m) * 1.0 / p[k];
		}
	}
	// real actual matrix multiplication
	auto tmp = RtP * R;
	return (double)nu * tmp;
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
				std::vector sigma     = calculateSigmaFromDeterministicPhi(method[l]);
				std::vector<double> p = sigma;
				// now normalize by row sum!
				coretools::normalize(p);
				// add to JP
				fisherInfo += _calcFisherInfo(Gamma, p, sigma, method[l].nu(index_in_counts), method[l]);
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
				_curLL[ix_linNumMethLoc] =
				    _distr->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc, i);
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
				_tryLL[ix_linNumMethLoc] =
				    _distr->calculateLL(sigma, method[l], p, index_in_counts, ix_linNumMethLoc, i);
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
			_tryLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc, i);
		}
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
};

double TModelBase::calculateLLRatio_perMethod(size_t i, const TData &Data) {
	// loop over all locations and timepoints of a particular method of a single species
	const auto &method = Data[i];

	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { // chame gd hei ga
		return 0.0;
	}

	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		const auto &loc        = method[j];
		size_t index_in_counts = method.getIndexInCounts(_speciesIDinUniqueContainer);
		std::vector sigma      = calculateSigma(loc);
		std::vector<double> p  = sigma;
		// now normalize by row sum!
		coretools::normalize(p);
		// calculate Likelihood
		size_t ix_linNumMethLoc  = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
		_tryLL[ix_linNumMethLoc] = _distr->calculateLL(sigma, method[j], p, index_in_counts, ix_linNumMethLoc, i);
	}
	return coretools::containerSum(_tryLL) - coretools::containerSum(_curLL);
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
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
		std::swap(_tryLL[ix_linNumMethLoc], _curLL[ix_linNumMethLoc]);
	}
}

double TModelBase::curLL(size_t ij) const { return _curLL[ij]; }

double TModelBase::curLLPerMethod(size_t i, const TData &Data) {
	// used when a method has changed
	const auto &method = Data[i];
	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return 0.0; }
	double LL = 0.0;
	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
		LL += _curLL[ix_linNumMethLoc];
	}
	return LL;
}

double TModelBase::sumCurLL() const { return coretools::containerSum(_curLL); }

void TModelBase::setCurLL(size_t ij, double Value) { _curLL[ij] = Value; }

void TModelBase::setCurLLPerMethod(const std::vector<double> &Values, size_t i, const TData &Data) {
	// used when a method has changed
	const auto &method = Data[i];
	if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }
	for (size_t j = 0; j < method.size(); ++j) { // loop over locations
		size_t ix_linNumMethLoc = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
		_curLL[ix_linNumMethLoc]               = Values[j];
	}
}

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

void TDeterministic::estimateInitial(const TData &Data, const TUniqueContainer<size_t> &Timepoints) {
	_estimateInitialParameters(Data);
	_gamma.estimateInitialGamma(Data, Timepoints, _speciesIDinUniqueContainer);
	_fillPhi();
	_oldPhi = _newPhi;
}

void TDeterministic::initializeTempVariables(const TData &Data) {
	_gamma.initializeJP(Data, *this);
	_initializeLL(Data);
}

void TDeterministic::_fillPhi() { _newPhi = _gamma.phi(); }

TypePhi TDeterministic::phi(size_t /*location_id*/, size_t timepoint_id) const { return _newPhi[timepoint_id]; }

TypePhi TDeterministic::deterministicPhi(size_t timepoint_id) const { return _newPhi[timepoint_id]; }

void TDeterministic::update(const TData &Data) {
	// update distribution parameters
	_updateDistributionParameters(Data);

	for (size_t m = 0; m < _gamma.numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma.updateSpecificIndex(m)) {
			_oldPhi = _newPhi;
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

	// update gamma summaries
	_gamma.updateSummaryGammaPosterior();
}

void TDeterministic::simulate(TData &Data, const std::vector<double> &n_bar, bool SimulateCovariates) {

	_gamma.setSimulatedGammaToZero();
	_fillPhi();
	_simulateDistributionParameters(Data, n_bar, SimulateCovariates);
};

// write posterior summaries of gamma's
void TDeterministic::resetSummaryGammaPosterior() { _gamma.resetSummaryGammaPosterior(); }

void TDeterministic::writeSummaryGammaPosterior(const std::string &Prefix) {
	_gamma.writeSummaryGammaPosterior(Prefix);
}

//------------------------------------------
// TStochastic
//------------------------------------------

TStochastic::TStochastic(size_t SpeciesIDinUniqueContainer, const TData &Data, size_t NumLocations,
                         size_t NumTimePoints, const std::unique_ptr<TBirpNegBinomAddOn> &NegBinomAddOn,
                         stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                         stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                         stattools::TParameterTyped<TypeLogPhi, 2> *LogPhi,
                         stattools::TParameterTyped<TypeLogSigma, 1> *LogSigma)
    : TModelBase(SpeciesIDinUniqueContainer, Data, NegBinomAddOn), _prior(LogPhi, LogSigma, Gamma, TimesOfChange),
      _logPhi(LogPhi), _logSigma(LogSigma) {

	_jointUpdater = std::make_unique<stattools::TUpdateUniqueIrregular<TypeLogPhi>>(
	    NumLocations * NumTimePoints, true, stattools::ProposalKernel::normal, "jointLogPhi");
	_jointUpdater->setProposalWidth("0.1");
}

void TStochastic::initialize(const TData &Data, const TUniqueContainer<std::string> &Locations,
                             const TUniqueContainer<std::size_t> &TimePoints) {
	_initializeDistribution(Data, Locations);
	_prior.initialize(TimePoints);

	// create name classes for locations
	auto locNames = std::make_shared<coretools::TNamesStrings>(Locations.vec());
	std::vector<std::string> tmp(TimePoints.size());
	for (size_t t = 0; t < TimePoints.size(); ++t) { tmp[t] = coretools::str::toString(TimePoints[t]); }
	auto timeNames = std::make_shared<coretools::TNamesStrings>(tmp);

	_logPhi->initStorage({Locations.size(), TimePoints.size()}, {locNames, timeNames});
	_logSigma->initStorage();
}

void TStochastic::_setPhiFromDeterministic() {
	auto deterministicPhi = _prior.calculateDeterministicPhi();

	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 0; k < _logPhi->dimensions()[1]; ++k) { // loop over timepoints
			const auto index = _logPhi->getIndex({j, k});       // get linear index
			_logPhi->setInitial(index, log(deterministicPhi[k]));
		}
	}
}

void TStochastic::estimateInitial(const TData &Data, const TUniqueContainer<size_t> &Timepoints) {
	_estimateInitialParameters(Data);
	_prior.estimateInitialGamma(Data, Timepoints, _speciesIDinUniqueContainer);

	// initialize phi based on deterministic phi
	_setPhiFromDeterministic();

	// set sigma
	_logSigma->setInitial(0.0);
}

void TStochastic::initializeTempVariables(const TData &Data) {
	_prior.initializeTempVariables(Data, *this, true);
	_initializeLL(Data);
}

void TStochastic::_updateLogPhi(size_t j, size_t k, const TData &Data) {
	// propose a new value
	const size_t ix = _logPhi->getIndex({j, k});
	_logPhi->updateSpecificIndex(ix);

	const auto logPriorRatio = _prior.calculateLogPriorRatio(j, k) + _prior.calculateLogPriorRatio(j, k + 1);
	// only need to calculate LL for this location j and timepoint k, but across all methods (if they have data for j
	// and k)
	const auto LLRatio       = calculateLLRatio_perLocation(j, Data);
	const double logH        = LLRatio + logPriorRatio;
	if (_logPhi->acceptOrReject(logH, ix)) {
		_prior.swapTryCur(j, k);
		_prior.swapTryCur(j, k + 1);
		swapTryCur_perLocation(j, Data);
	}
}

void TStochastic::_updateLogPhi(const TData &Data) {
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) {
			// loop over all timepoints (ignore first, since this is always 1)
			_updateLogPhi(j, k, Data);
		}
	}
}

void TStochastic::_jointUpdateLogPhiGamma(const TData &Data) {
	if (!_prior.gammaIsUpdated()) { return; }
	for (size_t m = 0; m < _prior.numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		// propose new gamma
		const auto Delta = _prior.updateGamma_forLogPhiJointUpdate(m);
		// adjust logPhi by Delta_k
		for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) {     // loop over all locations
			for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) { // loop over all timepoints
				const size_t linear     = _logPhi->getIndex({j, k});
				const double logPhi_new = _logPhi->value(linear) + Delta[k];
				_logPhi->set(linear, logPhi_new);
			}
		}

		// calculate LL: no need to add prior on phi (= transition probabilities), since they didn't change!
		const double LLRatio = calculateLLRatio(Data);
		if (_prior.evaluateGamma_forLogPhiJointUpdate(m, LLRatio, Data, *this)) {
			// accepted: swap LL of all locations
			swapTryCur();
		} else {
			// rejected: re-set all logPhi
			for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) {     // loop over all locations
				for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) { // loop over all timepoints
					_logPhi->reset(_logPhi->getIndex({j, k}));
				}
			}
		}
	}
	// update summaries of gamma posterior after update has finished
	_prior.updateSummaryGammaPosterior();
}

void TStochastic::_shiftLogPhis(size_t j, size_t k, double diff) {
	for (size_t l = k; l < _logPhi->dimensions()[1]; ++l) {
		const size_t ix = _logPhi->getIndex({j, l});
		_logPhi->set(ix, _logPhi->value(ix) + diff);
	}
}

void TStochastic::_jointUpdateLogPhi(size_t j, size_t k, const TData &Data) {
	// propose a new value
	const size_t ix = _logPhi->getIndex({j, k});
	auto diff       = _jointUpdater->propose(0.0, ix);

	// calculate difference and add this difference to all logPhi on the right (including current phi)
	_shiftLogPhis(j, k, diff);

	// only transition probability from previous matters!
	const auto logPriorRatio = _prior.calculateLogPriorRatio(j, k);

	// only need to calculate LL for this location j and timepoint k, but across all methods (if they have data for j
	// and k)
	const auto LLRatio  = calculateLLRatio_perLocation(j, Data);
	const double logH   = LLRatio + logPriorRatio;
	const size_t ixLast = _logPhi->getIndex({j, _logPhi->dimensions()[1]});
	coretools::TRange range(ix, ixLast, 1);
	if (_logPhi->acceptOrReject<false>(logH, range)) {
		_prior.swapTryCur(j, k);
		swapTryCur_perLocation(j, Data);
	} else {
		_jointUpdater->reject(ix);
	}
}

void TStochastic::_jointUpdateLogPhi(const TData &Data) {
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) {
			// loop over all timepoints (ignore first, since this is always 1)
			_jointUpdateLogPhi(j, k, Data);
		}
	}
}

void TStochastic::_updateLogSigma() {
	// propose a new value
	if (!_logSigma->update()) { return; }

	double LLRatio = 0.0;
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) { // loop over all locations
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) { LLRatio += _prior.calculateLogPriorRatio(j, k); }
	}

	const double logH = LLRatio + _logSigma->getLogPriorRatio();
	if (_logSigma->acceptOrReject(logH, 0)) {
		// accepted: update transition probabilities
		_prior.swapTryCur();
	}
}

void TStochastic::update(const TData &Data) {
	// update distribution parameters
	_updateDistributionParameters(Data);

	// update gamma
	_prior.updateGamma(Data, *this);

	// update logPhi
	if (_logPhi->isUpdated()) {
		_updateLogPhi(Data);
		_jointUpdateLogPhiGamma(Data);
		_jointUpdateLogPhi(Data);
	}

	// update logSigma
	if (_logSigma->isUpdated()) { _updateLogSigma(); }
}

void TStochastic::simulate(TData &Data, const std::vector<double> &n_bar, bool SimulateCovariates) {

	// simulate gamma
	_prior.setSimulatedGammaToZero();
	// set logSigma to something small (only if it has not been specified by the user)
	_logSigma->setInitial(-1.0);

	_prior.initializeTempVariables(Data, *this, false);

	// simulate phi from transition probabilities
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) { // loop over all locations
		size_t ix = _logPhi->getIndex({j, 0});
		_logPhi->setInitial(ix, 0.0); // set first phi = 1 -> logPhi = 0
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) {
			ix = _logPhi->getIndex({j, k});
			_logPhi->setInitial(ix, _prior.sampleFromTransitionProbability(j, k));
		}
	}

	_simulateDistributionParameters(Data, n_bar, SimulateCovariates);
};

coretools::StrictlyPositive<double> TStochastic::phi(size_t j, size_t k) const {
	const double phi = exp(_logPhi->value(_logPhi->getIndex({j, k})));
	if (phi < coretools::StrictlyPositive<double>::min()) { return coretools::StrictlyPositive<double>::min(); }
	if (phi > coretools::StrictlyPositive<double>::max()) { return coretools::StrictlyPositive<double>::max(); }
	return phi;
}

TypePhi TStochastic::deterministicPhi(size_t timepoint_id) const { return _prior.deterministicPhi(timepoint_id); }

// write posterior summaries of gamma's
void TStochastic::resetSummaryGammaPosterior() { _prior.resetSummaryGammaPosterior(); }

void TStochastic::writeSummaryGammaPosterior(const std::string &Prefix) { _prior.writeSummaryGammaPosterior(Prefix); }
