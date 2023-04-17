//
// Created by caduffm on 4/13/23.
//

#include "TDistributions.h"
#include "TLikelihood.h"
#include "stattools/ParametersObservations/TIndexSampler.h"

using namespace coretools::instances;

//--------------------------------
// TPoissonEpochModel
//--------------------------------

TPoissonDistribution::TPoissonDistribution() : TDistributionBase() {}

double TPoissonDistribution::calculateLL(const std::vector<double> & /*sigma*/, const TLocations &loc,
                                         const std::vector<double> &p, size_t index_in_counts,
                                         size_t /*ix_linNumMethLoc*/) const {
	double LL = 0.0;
	for (size_t k = 0; k < loc.size(); k++) { LL += (double)loc[k].counts_per_species(index_in_counts) * log(p[k]); }
	return LL;
}

void TPoissonDistribution::simulate(
    TData &Data, const TModelBase &Model, size_t N0,
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>> &detectionProbs,
    const TUniqueContainer<std::string> &locations, bool SimulateCovariates) {
	for (auto &method : Data) { method.simulatePoisson(N0, Model, detectionProbs, locations, SimulateCovariates); }
}

//--------------------------------
// TNegBinEpochModel
//--------------------------------

TNegBinDistribution::TNegBinDistribution(stattools::TParameterTyped<TypeMu, 1> *Mu,
                                         stattools::TParameterTyped<TypeLogA, 1> *LogA)
    : TDistributionBase() {
	_mu   = Mu;
	_logA = LogA;

	// set max of logA to something reasonable: otherwise we have numeric issues with gammaLog function
	TypeLogA::setMax(10.0);
}

void TNegBinDistribution::initializeInferred(const TUniqueContainer<size_t> &Timepoints,
                                             const TUniqueContainer<std::string> &Locations, const TData &Data) {
	// initialize mu
	auto muNames = std::make_shared<coretools::TNamesStrings>();
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &m = Data[i];
		if (m.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			for (size_t j = 0; j < m.size(); ++j) {
				const std::string name = m.name() + "_" + Locations[m[j].location_id()];
				muNames->addName({name});
			}
		}
	}
	_mu->initStorage({Data.numMethLoc(_speciesIDinUniqueContainer)}, {muNames});

	// initialize a: single element
	_logA->initStorage({1});
}

double TNegBinDistribution::calculateLL(const std::vector<double> &sigma, const TLocations &loc,
                                        const std::vector<double> & /*p*/, size_t index_in_counts,
                                        size_t ix_linNumMethLoc) const {
	double sum_a   = 0.0;
	double term2   = 0.0;
	TypeNu n       = loc.nu(index_in_counts);
	const double a = exp(_logA->value());
	for (size_t k = 0; k < loc.size(); k++) {  // loop over timepoints
		double alpha_k = a * _mu->value(ix_linNumMethLoc) * sigma[k];
		sum_a += alpha_k;                      // summing up alphas
		term2 += coretools::gammaLog((double)loc[k].counts_per_species(index_in_counts) + alpha_k) -
		         coretools::gammaLog(alpha_k); // x_k + 1 drops out in hastings ratio
	}
	return coretools::gammaLog(sum_a) - coretools::gammaLog((double)n + sum_a) + term2; // (n+1) drops out
}

void TNegBinDistribution::_updateMu(const TData &Data, const TModelBase &Model) {
	stattools::TIndexSampler sampler(
	    Data.numMethLoc(_speciesIDinUniqueContainer)); // draw random indices from range [0, _numMethLoc]
	sampler.sampleIndices();

	for (size_t c = 0; c < sampler.length(); ++c) {
		// pick two indices i to update
		auto [c1, c2] = sampler.getIndexPair(c);

		auto [m_ix_c1, l_ix_c1] = Data.linear_to_i_j(_speciesIDinUniqueContainer, c1);
		auto [m_ix_c2, l_ix_c2] = Data.linear_to_i_j(_speciesIDinUniqueContainer, c2);
		const auto &m1          = Data[m_ix_c1];
		const auto &m2          = Data[m_ix_c2];
		const auto &loc1        = m1[l_ix_c1];
		const auto &loc2        = m2[l_ix_c2];

		_mu->updatePairWiseSumOneVector(c1, c2);

		std::vector<double> sigma1 = TLikelihood::calculateSigma(loc1, Model);
		std::vector<double> sigma2 = TLikelihood::calculateSigma(loc2, Model);

		size_t index_in_counts1 = m1.getIndexInCounts(_speciesIDinUniqueContainer);
		const auto newLL1       = calculateLL(sigma1, loc1, {}, index_in_counts1, c1);
		size_t index_in_counts2 = m2.getIndexInCounts(_speciesIDinUniqueContainer);
		const auto newLL2       = calculateLL(sigma2, loc2, {}, index_in_counts2, c2);

		const double LLRatio = newLL1 + newLL2 - _oldLL[c1] - _oldLL[c2]; // in log
		const double logH    = LLRatio + _mu->getLogPriorRatio(c1) + _mu->getLogPriorRatio(c2);

		coretools::TRange range(c1, c2);
		if (_mu->acceptOrReject(logH, range, c1)) {
			// accepted -> update LL
			_oldLL[c1] = newLL1;
			_oldLL[c2] = newLL2;
		}
	}
}

void TNegBinDistribution::_updateLogA(const TData &Data, const TModelBase &Model) {
	if (_logA->update()) {
		std::vector<double> newLL(Data.numMethLoc(_speciesIDinUniqueContainer));
		for (size_t ij = 0; ij < Data.numMethLoc(_speciesIDinUniqueContainer); ++ij) {
			auto [m_ix, l_ix]      = Data.linear_to_i_j(_speciesIDinUniqueContainer, ij);
			const auto &m          = Data[m_ix];
			const auto &loc        = m[l_ix];
			std::vector sigma      = TLikelihood::calculateSigma(loc, Model);
			size_t index_in_counts = m.getIndexInCounts(_speciesIDinUniqueContainer);
			newLL[ij]              = calculateLL(sigma, loc, {}, index_in_counts, ij);
		}
		const double LLRatio = coretools::containerSum(newLL) - coretools::containerSum(_oldLL);
		const double logH    = LLRatio + _logA->getLogPriorRatio();
		if (_logA->acceptOrReject(logH, 0)) {
			// accepted -> update LL
			_oldLL = newLL;
		}
	}
}

void TNegBinDistribution::update(const TData &Data, const TModelBase &Model) {
	if (_mu->isUpdated()) { _updateMu(Data, Model); }
	if (_logA->isUpdated()) { _updateLogA(Data, Model); }
}

void TNegBinDistribution::_estimateInitialMu(const TData &Data) {
	// initialize mu: relative nu's
	std::vector<double> nus(Data.numMethLoc(_speciesIDinUniqueContainer));
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) {
		auto [m_id, l_id] = Data.linear_to_i_j(_speciesIDinUniqueContainer, c);
		nus[c]            = std::max(coretools::UnsignedInt<uint32_t>(1),
		                             Data[m_id][l_id].nu(Data[m_id].getIndexInCounts(
                              _speciesIDinUniqueContainer))); // do not allow 0: problem for Dirichlet
	}
	coretools::normalize(nus);

	// initialize to very small value: not all i-j-combinations exists
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) { _mu->setInitial(c, nus[c]); }
	_mu->normalizeToSum1AfterUpdate();
}

void TNegBinDistribution::estimateInitialPriorParameters(const TData &Data) {
	// initialize mu
	_estimateInitialMu(Data);
	_logA->setInitial(0.0); // log(1) = 0
}

double TNegBinDistribution::_setAForSimulation(const TData &Data) {
	_logA->setInitial(log(Data.numMethLoc(_speciesIDinUniqueContainer)));
	// calculate a = exp(log(a))
	const double a = exp(_logA->value());
	return a;
}

void TNegBinDistribution::simulate(
    TData &Data, const TModelBase &Model, size_t N0,
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<coretools::Probability>> &detectionProbs,
    const TUniqueContainer<std::string> &locations, bool SimulateCovariates) {
	// set a
	double A = _setAForSimulation(Data);

	for (auto &method : Data) {
		if (!method.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }

		// simulate mus based on detectionProb and a: will return a vector with length = numLocations
		auto mus = method.simulateNB(N0, Model, detectionProbs, locations, A, SimulateCovariates);
		coretools::normalize(mus);
		for (size_t c = 0; c < method.size(); ++c) { _mu->set(c, mus[c]); }
	}
}
