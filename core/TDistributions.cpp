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

TNegBinDistribution::TNegBinDistribution(size_t SpeciesIDinUniqueContainer,
                                         stattools::TParameterTyped<TypeMuOrN, 1> *MuOrN,
                                         stattools::TParameterTyped<TypeA, 1> *A)
    : TDistributionBase(SpeciesIDinUniqueContainer) {
	_muOrN = MuOrN;
	_a     = A;

	_isMu = _muOrN->name() == "mu";
}

void TNegBinDistribution::_initializeMu(const TData &Data, const TUniqueContainer<std::string> &Locations) {
	auto names = std::make_shared<coretools::TNamesStrings>();
	for (size_t i = 0; i < Data.size(); ++i) {
		const auto &m = Data[i];
		if (m.hasDataForSpeciesID(_speciesIDinUniqueContainer)) {
			for (size_t j = 0; j < m.size(); ++j) {
				const std::string name = m.name() + "_" + Locations[m[j].location_id()];
				names->addName({name});
			}
		}
	}
	_muOrN->initStorage({Data.numMethLoc(_speciesIDinUniqueContainer)}, {names});
}

void TNegBinDistribution::_initializeN(const TUniqueContainer<std::string> &Locations) {
	auto names = std::make_shared<coretools::TNamesStrings>(Locations.vec());
	_muOrN->initStorage({Locations.size()}, {names});
}

void TNegBinDistribution::initialize(const TData &Data, const TUniqueContainer<std::string> &Locations) {
	// initialize mu
	if (_isMu) {
		_initializeMu(Data, Locations);
	} else { // N
		_initializeN(Locations);
	}

	// initialize a: number of methods
	auto names = std::make_shared<coretools::TNamesStrings>();
	for (size_t i = 0; i < Data.size(); ++i) { names->addName({Data[i].name()}); }
	_a->initStorage({Data.size()}, {names});
}

double TNegBinDistribution::_getMuOrN(size_t location_id, size_t ij) const {
	if (_isMu) {
		return _muOrN->value(ij); // if mu: method-location specific
	}
	return _muOrN->value(location_id); // if N: location-specific
}

double TNegBinDistribution::calculateLL(const std::vector<double> &sigma, const TLocations &loc,
                                        const std::vector<double> & /*p*/, size_t index_in_counts,
                                        size_t ix_linNumMethLoc, size_t ix_method) const {
	double alpha_0 = 0.0;
	double term2   = 0.0;
	const double a = _a->value(ix_method);
	for (size_t k = 0; k < loc.size(); k++) { // loop over timepoints
		const double mu = _getMuOrN(loc.location_id(), ix_linNumMethLoc);
		const double phi_s = sigma[k];
		double alpha_k = mu * phi_s / a;
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

			_muOrN->updatePairWiseSumOneVector(linear1, linear2);

			std::vector<double> sigma1 = Model.calculateSigma(loc1);
			std::vector<double> sigma2 = Model.calculateSigma(loc2);

			size_t index_in_counts = m.getIndexInCounts(_speciesIDinUniqueContainer);
			const auto newLL1      = calculateLL(sigma1, loc1, {}, index_in_counts, linear1, i);
			const auto newLL2      = calculateLL(sigma2, loc2, {}, index_in_counts, linear2, i);

			const double LLRatio = newLL1 + newLL2 - Model.curLL(linear1) - Model.curLL(linear2); // in log
			const double logH    = LLRatio + _muOrN->getLogPriorRatio(linear1) + _muOrN->getLogPriorRatio(linear2);

			coretools::TRange range(linear1, linear2);
			if (_muOrN->acceptOrReject(logH, range, linear1)) {
				// accepted -> update LL
				Model.setCurLL(linear1, newLL1);
				Model.setCurLL(linear2, newLL2);
			}
		}
	}
}

std::pair<double, std::vector<double>> TNegBinDistribution::_calcHastingsRatioN(const std::vector<size_t> &ij,
                                                                                const TData &Data, TModelBase &Model) {
	std::vector<double> newLL(ij.size());
	double LLRatio = 0.0;
	for (size_t c = 0; c < ij.size(); ++c) {
		// loop over methods-locations combinations relevant for second location
		const auto [i, j]         = Data.linear_to_i_j(_speciesIDinUniqueContainer, ij[c]);
		size_t index_in_counts    = Data[i].getIndexInCounts(_speciesIDinUniqueContainer);
		const auto &loc           = Data[i][j];
		std::vector<double> sigma = Model.calculateSigma(loc);
		newLL[c]                  = calculateLL(sigma, loc, {}, index_in_counts, ij[c], i);
		LLRatio += newLL[c] - Model.curLL(ij[c]);
	}
	return {LLRatio, newLL};
}

void TNegBinDistribution::_updateN(const TData &Data, TModelBase &Model) {
	// sum_j N_j = 1
	// draw random indices from range [0, J-1] (numLocations)
	stattools::TIndexSampler sampler(_muOrN->totalSize());
	sampler.sampleIndices();

	for (size_t c = 0; c < sampler.length(); ++c) {
		// pick two indices j to update: correspond to indices in unique container
		auto [location_id_1, location_id_2] = sampler.getIndexPair(c);
		_muOrN->updatePairWiseSumOneVector(location_id_1, location_id_2);

		// for those location-ids: which ones are the relevant linear indices (ij)?
		// -> all the methods where this location occurs
		const auto &ij_1 = Data.locationId_to_ij(_speciesIDinUniqueContainer, location_id_1);
		const auto &ij_2 = Data.locationId_to_ij(_speciesIDinUniqueContainer, location_id_2);

		const auto [LLRatio_1, newLL_1] = _calcHastingsRatioN(ij_1, Data, Model);
		const auto [LLRatio_2, newLL_2] = _calcHastingsRatioN(ij_2, Data, Model);

		const double logH =
		    LLRatio_1 + LLRatio_2 + _muOrN->getLogPriorRatio(location_id_1) + _muOrN->getLogPriorRatio(location_id_2);
		coretools::TRange range(location_id_1, location_id_2);
		if (_muOrN->acceptOrReject(logH, range, location_id_1)) {
			// accepted -> update LL
			for (size_t a = 0; a < ij_1.size(); ++a) { Model.setCurLL(ij_1[a], newLL_1[a]); }
			for (size_t a = 0; a < ij_2.size(); ++a) { Model.setCurLL(ij_2[a], newLL_2[a]); }
		}
	}
}

void TNegBinDistribution::_updateA(const TData &Data, TModelBase &Model) {
	for (size_t i = 0; i < Data.size(); ++i) { // loop over methods
		if (!Data[i].hasDataForSpeciesID(_speciesIDinUniqueContainer)) { continue; }
		_a->updateSpecificIndex(i);
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
		const double logH    = LLRatio + _a->getLogPriorRatio(i);
		if (_a->acceptOrReject(logH, i)) {
			// accepted -> update LL
			Model.setCurLLPerMethod(newLL, i, Data);
		}
	}
}

void TNegBinDistribution::update(const TData &Data, TModelBase &Model) {
	if (_muOrN->isUpdated()) {
		if (_isMu) {
			_updateMu(Data, Model);
		} else {
			_updateN(Data, Model);
		}
	}
	if (_a->isUpdated()) { _updateA(Data, Model); }
}

std::vector<double> TNegBinDistribution::_getNuPerIJ(const TData &Data) {
	std::vector<double> nus(Data.numMethLoc(_speciesIDinUniqueContainer));
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) {
		auto [m_id, l_id] = Data.linear_to_i_j(_speciesIDinUniqueContainer, c);
		nus[c]            = std::max(coretools::UnsignedInt<uint64_t>(1),
		                             Data[m_id][l_id].nu(Data[m_id].getIndexInCounts(
                              _speciesIDinUniqueContainer))); // do not allow 0: problem for Dirichlet
	}
	return nus;
}

void TNegBinDistribution::_estimateInitialMu(const TData &Data) {
	// initialize mu: relative nu's
	auto nu = _getNuPerIJ(Data);
	for (size_t c = 0; c < Data.numMethLoc(_speciesIDinUniqueContainer); ++c) { _muOrN->setInitial(c, nu[c]); }

	// now normalize mu
	for (size_t i = 0; i < Data.size(); ++i) { // loop over methods
		double sum = 0.0;
		for (size_t j = 0; j < Data[i].size(); ++j) { // loop over locations and sum up values of mu
			const double mu = _muOrN->value(Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j));
			sum += mu;
		}
		for (size_t j = 0; j < Data[i].size(); ++j) { // loop over locations and divide by sum
			const size_t ix = Data.i_j_to_linear(_speciesIDinUniqueContainer, i, j);
			const double mu = _muOrN->value(ix);
			_muOrN->setInitial(ix, mu / sum);
		}
	}
}

void TNegBinDistribution::_estimateInitialN(const TData &Data) {
	// initialize N: relative nu's
	auto nu = _getNuPerIJ(Data);

	std::vector<double> nu_per_location(_muOrN->totalSize(), 0.0);
	for (size_t location_id = 0; location_id < _muOrN->totalSize(); ++location_id) {
		const auto &ij = Data.locationId_to_ij(_speciesIDinUniqueContainer, location_id);
		for (auto c : ij) { nu_per_location[location_id] += nu[c]; }
	}

	// normalize
	coretools::normalize(nu_per_location);
	for (size_t location_id = 0; location_id < _muOrN->totalSize(); ++location_id) {
		_muOrN->setInitial(location_id, nu_per_location[location_id]);
	}
}

void TNegBinDistribution::estimateInitialParameters(const TData &Data) {
	// initialize mu
	if (_isMu) {
		_estimateInitialMu(Data);
	} else {
		_estimateInitialN(Data);
	}
	for (size_t i = 0; i < Data.size(); ++i) { _a->setInitial(i, 1.0); }
}

std::vector<double> TNegBinDistribution::_setAForSimulation(const TData &Data) {
	std::vector<double> a(Data.size());
	for (size_t i = 0; i < Data.size(); ++i) {
		_a->setInitial(i, (double)Data[i].size());
		a[i] = _a->value(i);
	}
	return a;
}

void TNegBinDistribution::simulate(TData &Data, const TModelBase &Model,
                                   const std::vector<double> &AvgNumObsPerTimepoint, bool SimulateCovariates) {
	// set a
	auto A = _setAForSimulation(Data);

	for (size_t i = 0; i < Data.size(); ++i) {
		if (!Data[i].hasDataForSpeciesID(_speciesIDinUniqueContainer)) { return; }

		Data[i].simulateNB(AvgNumObsPerTimepoint[i], Model, A[i], SimulateCovariates);

		// set mu: 1 / J_i (all locations contribute equally)
		if (_isMu) {
			for (size_t c = 0; c < Data[i].size(); ++c) {
				_muOrN->set(Data.i_j_to_linear(_speciesIDinUniqueContainer, i, c), 1.0 / (double)Data[i].size());
			}
		}
	}
	// set N: 1 / J (all locations contribute equally)
	if (!_isMu) {
		for (size_t j = 0; j < _muOrN->totalSize(); ++j) { _muOrN->set(j, 1.0 / (double)_muOrN->totalSize()); }
	}
}
