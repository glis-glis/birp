//
// Created by caduffm on 4/13/23.
//

#include "TDistributions.h"
#include "TModels.h"
#include "coretools/devtools.h"
#include "stattools/ParametersObservations/TIndexSampler.h"

using namespace coretools::instances;

//--------------------------------
// TDistributionBase
//--------------------------------

TDistributionBase::TDistributionBase(size_t SpeciesIDinUniqueContainer)
    : _speciesIDinUniqueContainer(SpeciesIDinUniqueContainer) {}

//--------------------------------
// TPoissonDistribution
//--------------------------------

TPoissonDistribution::TPoissonDistribution(size_t SpeciesIDinUniqueContainer)
    : TDistributionBase(SpeciesIDinUniqueContainer) {}

double TPoissonDistribution::calculateLL(const std::vector<double> & /*sigma*/, const TLocations &loc,
                                         const std::vector<double> &p, size_t index_in_counts,
                                         size_t /*ix_linNumMethLoc*/, size_t /*ix_method*/) const {
	double LL = 0.0;
	for (size_t k = 0; k < loc.size(); k++) { LL += (double)loc[k].counts_per_species(index_in_counts) * log(p[k]); }
	return LL;
}

void TPoissonDistribution::simulate(TData &Data, const TModelBase &Model,
                                    const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) {

	for (size_t i = 0; i < Data.size(); ++i) {
		if (!Data[i].hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }
		Data[i].simulatePoisson(AvgNumObsPerTimepoint[i], Model, SimulateCovariates);
	}
}

//--------------------------------
// TNegBinDistribution
//--------------------------------

TNegBinDistribution::TNegBinDistribution(size_t SpeciesIDinUniqueContainer, stattools::TParameterTyped<TypeMu, 1> *Mu,
                                         stattools::TParameterTyped<TypeLogA, 1> *LogA)
    : TDistributionBase(SpeciesIDinUniqueContainer) {
	_mu   = Mu;
	_logA = LogA;

	// set max of logA to something reasonable: otherwise we have numeric issues with gammaLog function
	TypeLogA::setMax(10.0);
}

void TNegBinDistribution::initialize(const TData &Data, const TUniqueContainer<std::string> &Locations) {
	// initialize mu
	auto muNames     = std::make_shared<coretools::TNamesStrings>();
	auto methodNames = std::make_shared<coretools::TNamesStrings>();
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &m = Data[i];
		methodNames->addName({m.name()});
		if (m.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			for (size_t j = 0; j < m.size(); ++j) {
				const std::string name = m.name() + "_" + Locations[m[j].location_id()];
				muNames->addName({name});
			}
		}
	}
	_mu->initStorage({Data.numMethLoc(_speciesIDinUniqueContainer)}, {muNames});

	// initialize a: number of methods
	_logA->initStorage({Data.size()}, {methodNames});
}

double TNegBinDistribution::calculateLL(const std::vector<double> &sigma, const TLocations &loc,
                                        const std::vector<double> & /*p*/, size_t index_in_counts,
                                        size_t ix_linNumMethLoc, size_t ix_method) const {
	double alpha_0 = 0.0;
	double term2   = 0.0;
	const double a = exp(_logA->value(ix_method));
	for (size_t k = 0; k < loc.size(); k++) { // loop over timepoints
		double alpha_k = a * _mu->value(ix_linNumMethLoc) * sigma[k];
		alpha_0 += alpha_k;
		const double x_k = (double)loc[k].counts_per_species(index_in_counts);
		term2 += coretools::diffGammaLog(alpha_k, x_k);
	}
	const double nu = (double)loc.nu(index_in_counts);
	return coretools::diffGammaLog(alpha_0, nu) - term2;
}

void TNegBinDistribution::_updateMu(const TData &Data, TModelBase &Model) {
	// sum_j mu_ij = 1
	for (size_t i = 0; i < Data.size(); ++i) { // loop over methods
		// if no data for this species: skip method
		const auto &m = Data[i];
		if (!m.hasDataForSpeciesID(_speciesIDinUniqueContainer)) { continue; }
		// draw random indices from range [0, J_i] (numLocations for current method)
		stattools::TIndexSampler sampler(m.size());
		sampler.sampleIndices();

		for (size_t c = 0; c < sampler.length(); ++c) {
			// pick two indices j to update
			auto [j1, j2] = sampler.getIndexPair(c);

			const size_t linear1 = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j1);
			const size_t linear2 = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j2);
			const auto &loc1     = m[j1];
			const auto &loc2     = m[j2];

			_mu->updatePairWiseSumOneVector(linear1, linear2);

			std::vector<double> sigma1 = Model.calculateSigma(loc1);
			std::vector<double> sigma2 = Model.calculateSigma(loc2);

			size_t index_in_counts = m.getIndexInCounts(_speciesIDinUniqueContainer);
			const auto newLL1      = calculateLL(sigma1, loc1, {}, index_in_counts, linear1, i);
			const auto newLL2      = calculateLL(sigma2, loc2, {}, index_in_counts, linear2, i);

			const double LLRatio = newLL1 + newLL2 - Model.curLL(linear1) - Model.curLL(linear2); // in log
			const double logH    = LLRatio + _mu->getLogPriorRatio(linear1) + _mu->getLogPriorRatio(linear2);

			coretools::TRange range(linear1, linear2);
			if (_mu->acceptOrReject(logH, range, linear1)) {
				// accepted -> update LL
				Model.setCurLL(linear1, newLL1);
				Model.setCurLL(linear2, newLL2);
			}
		}
	}
}

void TNegBinDistribution::_updateLogA(const TData &Data, TModelBase &Model) {
	for (size_t i = 0; i < Data.size(); ++i) { // loop over methods
		if (!Data[i].hasDataForSpeciesID(_speciesIDinUniqueContainer)) { continue; }
		_logA->updateSpecificIndex(i);
		const auto &m = Data[i];
		std::vector<double> newLL(m.size());
		for (size_t j = 0; j < m.size(); ++j) {
			const auto &loc              = m[j];
			std::vector sigma            = Model.calculateSigma(loc);
			const size_t index_in_counts = m.getIndexInCounts(_speciesIDinUniqueContainer);
			const size_t ij              = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			newLL[j]                     = calculateLL(sigma, loc, {}, index_in_counts, ij, i);
		}
		const double LLRatio = coretools::containerSum(newLL) - Model.curLLPerMethod(i, Data);
		const double logH    = LLRatio + _logA->getLogPriorRatio(i);
		if (_logA->acceptOrReject(logH, i)) {
			// accepted -> update LL
			Model.setCurLLPerMethod(newLL, i, Data);
		}
	}
}

void TNegBinDistribution::update(const TData &Data, TModelBase &Model) {
	if (_mu->isUpdated()) { _updateMu(Data, Model); }
	if (_logA->isUpdated()) { _updateLogA(Data, Model); }
}

void TNegBinDistribution::_estimateInitialMu(const TData &Data) {
	// initialize mu: relative nu's
	std::vector<double> nus(Data.numMethLoc(_speciesIDinUniqueContainer));
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) {
		auto [m_id, l_id] = Data.linear_to_i_j(_speciesIDinUniqueContainer, c);
		nus[c]            = std::max(coretools::UnsignedInt<uint64_t>(1),
		                             Data[m_id][l_id].nu(Data[m_id].getIndexInCounts(
                              _speciesIDinUniqueContainer))); // do not allow 0: problem for Dirichlet
	}
	coretools::normalize(nus);

	// initialize to very small value: not all i-j-combinations exists
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) { _mu->setInitial(c, nus[c]); }
	// now normalize mu
	for (size_t i = 0; i < Data.size(); ++i) { // loop over methods
		double sum = 0.0;
		for (size_t j = 0; j < Data[i].size(); ++j) { // loop over locations and sum up values of mu
			const double mu = _mu->value(Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j));
			sum += mu;
		}
		for (size_t j = 0; j < Data[i].size(); ++j) { // loop over locations and divide by sum
			const size_t ix = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			const double mu = _mu->value(ix);
			_mu->setInitial(ix, mu / sum);
		}
	}
}

void TNegBinDistribution::estimateInitialParameters(const TData &Data) {
	// initialize mu
	_estimateInitialMu(Data);
	for (size_t i = 0; i < Data.size(); ++i) { _logA->setInitial(i, 0.0); } // log(1) = 0
}

std::vector<double> TNegBinDistribution::_setAForSimulation(const TData &Data) {
	std::vector<double> a(Data.size());
	for (size_t i = 0; i < Data.size(); ++i) {
		_logA->setInitial(i, log(Data[i].size()));
		a[i] = exp(_logA->value(i)); // calculate a = exp(log(a))
	}
	return a;
}

void TNegBinDistribution::simulate(TData &Data, const TModelBase &Model,
                                   const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) {
	// set a
	auto A = _setAForSimulation(Data);

	for (size_t i = 0; i < Data.size(); ++i) {
		if (!Data[i].hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }

		// simulate mus based on detectionProb and a: will return a vector with length = numLocations
		auto mus = Data[i].simulateNB(AvgNumObsPerTimepoint[i], Model, A[i], SimulateCovariates);
		coretools::normalize(mus);
		for (size_t c = 0; c < Data[i].size(); ++c) { _mu->set(c, mus[c]); }
	}
}
