//
// Created by caduffm on 4/17/23.
//

#include "TGamma.h"

#include "TData.h"
#include "TModels.h"
#include "coretools/Distributions/TNormalDistr.h"
#include "coretools/Main/TParameters.h"
#include "coretools/devtools.h"

using namespace coretools::instances;

//------------------------------------------
// TOLSGamma
//------------------------------------------

TOLSGamma::TOLSGamma(size_t NumEpochs, const std::vector<TypeTime> &EpochStarts)
    : _numEpochs(NumEpochs), _epochStarts(EpochStarts) {
	_gamma.resize(_numEpochs, 0.0);
}

bool TOLSGamma::_timePointIsWithinEpoch(const TUniqueContainer<TypeTime> &Timepoints, size_t t, size_t m) const {
	// returns true if a timepoint is withing an epoch (borders included)
	return Timepoints[t] >= _epochStarts[m] && Timepoints[t] <= _epochStarts[m + 1];
}

auto TOLSGamma::_get_Y_m(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t m,
                         size_t SpeciesID) const {
	// calculates response variable y_m for each epoch as well as other variables needed for OLS
	std::vector<double> y_m;             // response variable
	std::vector<TypeTime> t_m;           // timepoints, same size as y_m
	std::vector<size_t> dataSeriesIndex; // index indicating which method-location-combi was used, same size as y_m
	std::vector<std::pair<size_t, size_t>> l_to_ij;  // map that indicates i-j for each linear l
	std::vector<size_t> numDataPointsWithZeroCounts; // count number of counts=0 for each linear l

	size_t L = 0;
	for (size_t i = 0; i < Data.size(); ++i) {        // loop over methods
		for (size_t j = 0; j < Data[i].size(); ++j) { // loop over locations
			bool dataSeries_has_Timepoint_within_epoch = false;
			size_t numDataPointsWithZeroCounts_this_ij = 0;
			for (size_t k = 0; k < Data[i][j].size(); ++k) { // loop over timepoints
				const auto &timepoint = Data[i][j][k];
				if (_timePointIsWithinEpoch(Timepoints, timepoint.timepoint_id(), m)) {
					// ignore datapoints with effort zero
					if (timepoint.effort() == 0.0) { continue; }

					dataSeries_has_Timepoint_within_epoch = true;

					// prevent issues with zero counts: log(0) = inf
					// exclude them from analysis, but remember how many were excluded
					size_t counts = timepoint.counts_per_species(SpeciesID);
					if (counts == 0) {
						++numDataPointsWithZeroCounts_this_ij;
						continue;
					}

					// calculate y
					const double y = log((double)counts / (double)timepoint.effort());

					// add to containers
					y_m.push_back(y);
					t_m.emplace_back(Timepoints[timepoint.timepoint_id()]);
					dataSeriesIndex.emplace_back(L);
				}
			}
			// check if method-location combination l had any data that falls within epoch
			// if all counts are zero for this combination l: can not do OLS, ignore
			// -> else: add entry to map l_to_ij that indicates i-j for a specific l
			if (dataSeries_has_Timepoint_within_epoch && numDataPointsWithZeroCounts_this_ij != Data[i][j].size()) {
				++L;
				l_to_ij.emplace_back(i, j);
				numDataPointsWithZeroCounts.emplace_back(numDataPointsWithZeroCounts_this_ij);
			}
		}
	}
	return std::make_tuple(y_m, t_m, dataSeriesIndex, l_to_ij, numDataPointsWithZeroCounts);
}

arma::vec TOLSGamma::_doOLS(const std::vector<double> &Y_m, const std::vector<TypeTime> &T_m,
                            const std::vector<size_t> &DataSeriesIndex, size_t L) {
	// fill matrix X
	arma::mat X(Y_m.size(), L + 1, arma::fill::zeros);
	for (size_t row = 0; row < X.n_rows; ++row) {
		const size_t l       = DataSeriesIndex[row];
		X(row, l)            = 1;
		X(row, X.n_cols - 1) = T_m[row]; // last column contains timepoints
	}

	// do OLS
	const arma::vec arma_y(Y_m);
	auto beta = arma::solve(X, arma_y);
	return beta.as_col(); // matrix -> vector
}

void TOLSGamma::_augmentDataSeriesIndexWithOLSEstimates(std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
                                                        size_t l, std::vector<std::vector<size_t>> &dataSeriesIndex,
                                                        size_t curEpoch, size_t otherEpoch) {
	// get method (i) and location (j) index of data series (l) of current epoch
	const size_t i                  = l_to_ij[curEpoch][l].first;
	const size_t j                  = l_to_ij[curEpoch][l].second;
	// check if other epoch (can be previous or next) already contains this specific combination of i and j
	const auto &l_to_ij_other_epoch = l_to_ij[otherEpoch];
	bool found                      = false;
	for (size_t l_other = 0; l_other < l_to_ij_other_epoch.size(); ++l_other) {
		if (l_to_ij_other_epoch[l_other].first == i && l_to_ij_other_epoch[l_other].second == j) {
			// method and location exists in previous epoch as well
			// -> use that index l to augment dataSeriesIndex
			// -> no need to augment l_to_ij since this ij already exists
			dataSeriesIndex[otherEpoch].push_back(l_other);
			found = true;
			break;
		}
	}
	// other epoch does not yet contain this specific combination of i and j
	// -> add a new entry {i, j} to l_to_ij
	// -> use this new l to augment dataSeriesIndex
	if (!found) {
		l_to_ij[otherEpoch].emplace_back(i, j);
		const size_t l_other = l_to_ij[otherEpoch].size() - 1;
		dataSeriesIndex[otherEpoch].push_back(l_other);
	}
}

bool TOLSGamma::_dontAugmentDataSeries(const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij,
                                       const std::vector<std::vector<size_t>> &dataSeriesIndex,
                                       const std::vector<std::vector<TypeTime>> &t, size_t l, size_t curEpoch,
                                       size_t otherEpoch, TypeTime EpochBorder) {
	// check if other epoch already has data at method i and location j and epoch border (start or end)
	// get method (i) and location (j) index of data series (l) of current epoch
	const size_t i = l_to_ij[curEpoch][l].first;
	const size_t j = l_to_ij[curEpoch][l].second;
	// loop over all timepoints of other epoch
	for (size_t r = 0; r < t[otherEpoch].size(); ++r) {
		// check if other epoch contains this method (i) and location (j) and the timepoint corresponding to EpochBorder
		const size_t other_l   = dataSeriesIndex[otherEpoch][r];
		const size_t other_i   = l_to_ij[otherEpoch][other_l].first;
		const size_t other_j   = l_to_ij[otherEpoch][other_l].second;
		const TypeTime other_t = t[otherEpoch][r];
		if (other_i == i && other_j == j && other_t == EpochBorder) { return true; }
	}
	return false;
}

void TOLSGamma::_augmentYWithOLSEstimatesStart(size_t curEpoch, const arma::vec &Beta,
                                               std::vector<std::vector<double>> &y,
                                               std::vector<std::vector<TypeTime>> &t,
                                               std::vector<std::vector<size_t>> &dataSeriesIndex,
                                               std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const {
	// start
	const TypeTime epochStart = _epochStarts[curEpoch];
	if (curEpoch == 0) {
		// first epoch: no need to augment previous epoch since there is none
		return;
	}
	for (size_t l = 0; l < l_to_ij[curEpoch].size() - 1; ++l) { // loop over all dataseries
		if (_dontAugmentDataSeries(l_to_ij, dataSeriesIndex, t, l, curEpoch, curEpoch - 1, epochStart)) { continue; }
		const auto alpha_l = Beta(l);
		const auto y_start = alpha_l + (double)epochStart * _gamma[curEpoch];
		// epoch start does not correspond to an existing timepoint -> augment data with y_start
		y[curEpoch - 1].push_back(y_start);
		t[curEpoch - 1].emplace_back(epochStart);
		_augmentDataSeriesIndexWithOLSEstimates(l_to_ij, l, dataSeriesIndex, curEpoch, curEpoch - 1);
	}
}

void TOLSGamma::_augmentYWithOLSEstimatesEnd(size_t curEpoch, const arma::vec &Beta,
                                             std::vector<std::vector<double>> &y, std::vector<std::vector<TypeTime>> &t,
                                             std::vector<std::vector<size_t>> &dataSeriesIndex,
                                             std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const {
	// end
	const TypeTime epochEnd = _epochStarts[curEpoch + 1];
	if (curEpoch == _numEpochs - 1) {
		// last epoch: no need to augment next epoch since there is none
		return;
	}
	for (size_t l = 0; l < l_to_ij[curEpoch].size() - 1; ++l) { // loop over all dataseries
		if (_dontAugmentDataSeries(l_to_ij, dataSeriesIndex, t, l, curEpoch, curEpoch + 1, epochEnd)) { continue; }
		const auto alpha_l = Beta(l);
		const auto y_end   = alpha_l + (double)epochEnd * _gamma[curEpoch];
		// epoch start does not correspond to an existing timepoint -> augment data with y_start
		y[curEpoch + 1].push_back(y_end);
		t[curEpoch + 1].emplace_back(epochEnd);
		_augmentDataSeriesIndexWithOLSEstimates(l_to_ij, l, dataSeriesIndex, curEpoch, curEpoch + 1);
	}
}

bool TOLSGamma::_epochIsIdentifiable(size_t L, const std::vector<TypeTime> &t,
                                     const std::vector<size_t> &dataSeriesIndex) {
	for (size_t l = 0; l < L; ++l) { // loop over all possible data series
		// count the number of occurrences of that particular data series (method/location combination) in the current
		// epoch
		size_t counter_l = 0;
		for (size_t r = 0; r < t.size(); ++r) { // loop over all data of current epoch
			if (dataSeriesIndex[r] == l) { ++counter_l; }
		}
		if (counter_l >= 2) { return true; }
	}
	return false;
}

bool TOLSGamma::_epochIsNotIdentifiableBecauseOfZeroCounts(size_t L, const std::vector<TypeTime> &t,
                                                           const std::vector<size_t> &dataSeriesIndex,
                                                           const std::vector<size_t> &numDataPointsWithZeroCounts) {
	if (numDataPointsWithZeroCounts.empty()) { return false; } // no data points at all!

	for (size_t l = 0; l < L; ++l) { // loop over all possible data series
		// count the number of occurrences of that particular data series (method/location combination) in the current
		// epoch
		size_t counter_l = 0;
		for (size_t r = 0; r < t.size(); ++r) { // loop over all data of current epoch
			if (dataSeriesIndex[r] == l) { ++counter_l; }
		}
		// add the number of data points that have zero counts for that particular data series
		counter_l += numDataPointsWithZeroCounts[l];
		if (counter_l >= 2) { return true; } // together with zero counts, we would have enough data!
	}
	return false; // we still don't have enough data
}

size_t TOLSGamma::_getNextEpochInitialGamma(const std::vector<bool> &EpochInitialized,
                                            const std::vector<std::vector<double>> &y,
                                            const std::vector<std::vector<std::pair<size_t, size_t>>> &l_to_ij) const {
	// find epoch that contains the most data points per degree of freedom (= number of parameters = L + 1)
	size_t curEpoch = 0;
	double max      = std::numeric_limits<double>::lowest();
	for (size_t m = 0; m < _numEpochs; ++m) {
		const size_t numDataPoints      = y[m].size();
		const size_t df                 = l_to_ij[m].size() + 1;
		const double numDataPointsPerDF = (double)numDataPoints / (double)df;
		if (!EpochInitialized[m] && numDataPointsPerDF > max) {
			curEpoch = m;
			max      = numDataPointsPerDF;
		}
	}
	return curEpoch;
}

std::vector<double> TOLSGamma::estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints,
                                                    size_t SpeciesID) {
	logfile().startIndent("Estimating initial gamma with OLS:");
	// fill response variable y
	std::vector<std::vector<double>> y(_numEpochs);
	std::vector<std::vector<TypeTime>> t(_numEpochs);
	std::vector<std::vector<size_t>> dataSeriesIndex(_numEpochs);
	std::vector<std::vector<std::pair<size_t, size_t>>> l_to_ij(_numEpochs);
	std::vector<std::vector<size_t>> numDataPointsWithZeroCounts(_numEpochs);
	for (size_t m = 0; m < _numEpochs; ++m) {
		std::tie(y[m], t[m], dataSeriesIndex[m], l_to_ij[m], numDataPointsWithZeroCounts[m]) =
		    _get_Y_m(Data, Timepoints, m, SpeciesID);
	}
	std::vector<bool> epochInitialized(_numEpochs, false);
	while (std::accumulate(epochInitialized.begin(), epochInitialized.end(), 0U) < _numEpochs) {
		// 1) find epoch that contains the most data points per degree of freedom (= number of parameters = L + 1)
		size_t curEpoch = _getNextEpochInitialGamma(epochInitialized, y, l_to_ij);
		logfile().listFlush("Epoch ", curEpoch + 1, " [", _epochStarts[curEpoch], ", ", _epochStarts[curEpoch + 1],
		                    "]:");

		// 2) we need at least one data series (method-location combination) with at least 2 timepoints
		const size_t L = l_to_ij[curEpoch].size();
		if (_epochIsIdentifiable(L, t[curEpoch], dataSeriesIndex[curEpoch])) {
			// 3) estimate intercept and gamma with OLS
			auto beta        = _doOLS(y[curEpoch], t[curEpoch], dataSeriesIndex[curEpoch], L);
			_gamma[curEpoch] = beta(L);

			// 4) calculate y_start and y_end and add it to existing data
			_augmentYWithOLSEstimatesStart(curEpoch, beta, y, t, dataSeriesIndex, l_to_ij);
			_augmentYWithOLSEstimatesEnd(curEpoch, beta, y, t, dataSeriesIndex, l_to_ij);
		} else {
			// check if epoch would be identifiable but isn't because we excluded counts = 0 for OLS
			if (_epochIsNotIdentifiableBecauseOfZeroCounts(L, t[curEpoch], dataSeriesIndex[curEpoch],
			                                               numDataPointsWithZeroCounts[curEpoch])) {
				// set gamma to zero: identifiable for MCMC, but not for OLS
				// don't augment y_start and y_end: not a reliable estimate of gamma
				_gamma[curEpoch] = 0.0;
			} else {
				UERROR("Need at least one method-location combination with two or more timepoints for epoch [",
				       _epochStarts[curEpoch], ", ", _epochStarts[curEpoch + 1], "]!");
			}
		}

		epochInitialized[curEpoch] = true;
		logfile().flush(" initial gamma = ", _gamma[curEpoch], "!");
		logfile().newLine();
	}
	logfile().endIndent();
	return _gamma;
}

//------------------------------------------
// TGamma
//------------------------------------------

TGamma::TGamma(TypeParamGamma *Gamma, TypeParamTOC *TimesOfChange) : _gamma(Gamma), _timesOfChange(TimesOfChange) {
	// find out how many epochs there are
	if (parameters().exists("timesOfChange")) {
		// Note: actual value of _timesOfChange will be set via stattools (name of parameter _timesOfChange is
		// "timesOfChange" -> parsed from command line)
		std::vector<std::string> tmp;
		parameters().fill("timesOfChange", tmp);
		logfile().list("Times of change = ", tmp, " (argument 'timesOfChange').");
		_numEpochs = tmp.size() + 1;
	} else {
		logfile().list("Will assume a single epoch (use argument 'timesOfChange' to change).");
		_numEpochs = 1;
	}
	// default: do not infer times of change
	if (!parameters().exists("inferTimesOfChange")) { _timesOfChange->getDefinition().update(false); }

	// calculate Jeffrey's prior?
	_noJP = parameters().exists("noJP");
	if (_noJP) { logfile().list("Will not use Jeffrey's prior on gamma (argument 'noJP')."); }
}

void TGamma::initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior) {
	_timesOfChange->initStorage(BirpPrior, {_numEpochs - 1}, {std::make_shared<coretools::TNamesIndices>()});
	_gamma->initStorage(BirpPrior, {_numEpochs}, {std::make_shared<coretools::TNamesIndices>()});

	_summaryGammaPosterior.resize(_numEpochs, std::vector<size_t>(_numEpochs, 0));

	_fillRho(Timepoints);
}

void TGamma::estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints, size_t SpeciesID) {
	if (_gamma->hasFixedInitialValue()) { return; }

	std::vector<TypeTime> epochStarts = getEpochStarts(Timepoints);
	TOLSGamma ols(_numEpochs, epochStarts);
	auto gamma = ols.estimateInitialGamma(Data, Timepoints, SpeciesID);
	assert(gamma.size() == _numEpochs);
	for (size_t m = 0; m < _numEpochs; ++m) { _gamma->set(m, gamma[m]); }
}

void TGamma::initializeJP(const TData &Data, const TModelBase &Model) {
	// set _oldLL and _oldJP
	_curJP = Model.calculateJeffreyPrior(Data, *this);
	_tryJP = _curJP;
}

std::vector<TypeTime> TGamma::getEpochStarts(const TUniqueContainer<TypeTime> &Timepoints) const {
	// return start of each epoch, including for a hypothetical epoch after the last epoch
	assert(!Timepoints.empty());

	auto minTime                      = *std::min(Timepoints.cbegin(), Timepoints.cend());
	std::vector<TypeTime> epochStarts = {minTime};
	if (_timesOfChange->size() != 0) {
		TypeTime maxTimeOfChange = _timesOfChange->max_value();
		TypeTime maxTime         = std::max((double)Timepoints.back(), maxTimeOfChange + 1.0);
		for (size_t i = 0; i < _timesOfChange->size(); i++) { epochStarts.push_back(_timesOfChange->value(i)); }
		epochStarts.push_back(maxTime);
	} else {
		epochStarts.push_back(*max_element(Timepoints.cbegin(), Timepoints.cend()));
	}
	return epochStarts;
}

void TGamma::_fillRho(const TUniqueContainer<TypeTime> &Timepoints) {
	std::vector<TypeTime> epochStarts = getEpochStarts(Timepoints);
	_rho.zeros((int)Timepoints.size(), (int)_numEpochs);

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

void TGamma::setSimulatedGammaToZero() {
	for (size_t g = 0; g < _gamma->size(); ++g) { _gamma->set(g, 0.0); }
}

// functions to update gamma

bool TGamma::updateSpecificIndex(size_t m) {
	if (_gamma->isUpdated()) {
		_gamma->propose(coretools::TRange(m));
		return true;
	}
	return false;
}

double TGamma::getLogPriorRatio(const TData &Data, const TModelBase &Model) {
	// -> prior on gamma: Jeffrey prior
	if (_noJP) { return 0.0; }
	_tryJP = Model.calculateJeffreyPrior(Data, *this);
	return _tryJP - _curJP;
}

bool TGamma::acceptOrReject(double logH, size_t m) {
	if (_gamma->acceptOrReject(logH, coretools::TRange(m))) {
		_curJP = _tryJP;
		return true;
	}
	return false;
}

void TGamma::resetSummaryGammaPosterior() {
	_counter = 0;
	for (size_t m = 0; m < _numEpochs; ++m) {
		std::fill(_summaryGammaPosterior[m].begin(), _summaryGammaPosterior[m].end(), 0);
	}
}

void TGamma::updateSummaryGammaPosterior() {
	++_counter;
	for (size_t m = 0; m < _numEpochs; ++m) {     // rows
		for (size_t n = 0; n < _numEpochs; ++n) { // column
			if (m == n) {
				// on the diagonal: count if gamma_m > 0
				if (_gamma->value(m) >= 0.0) { ++_summaryGammaPosterior[m][n]; }
			} else {
				// off the diagonal: count if gamma_m > gamma_n
				if (_gamma->value(m) > _gamma->value(n)) { ++_summaryGammaPosterior[m][n]; }
			}
		}
	}
}

void TGamma::writeSummaryGammaPosterior(const std::string &Prefix) {
	std::vector<std::string> rowNames(_numEpochs);
	for (size_t m = 0; m < _numEpochs; ++m) { rowNames[m] = "gamma_" + coretools::str::toString(m + 1); }
	auto header = rowNames;
	header.insert(header.begin(), "row_bigger_than_column");
	coretools::TOutputMaybeRcppFile file(Prefix + "_gammaSummaries.txt", header);
	for (size_t m = 0; m < _numEpochs; ++m) {
		file << rowNames[m];
		for (size_t n = 0; n < _numEpochs; ++n) { file << (double)_summaryGammaPosterior[m][n] / (double)_counter; }
		file.endln();
	}
}

// getters
size_t TGamma::numEpochs() const { return _numEpochs; }

size_t TGamma::numTimepoints() const { return _rho.n_rows; }

TypeGamma TGamma::gamma(size_t m) const { return _gamma->value(m); }

double TGamma::rho(size_t row, size_t col) const { return _rho(row, col); }

bool TGamma::gammaIsUpdated() const { return _gamma->isUpdated(); }

[[nodiscard]] std::vector<TypePhi> TGamma::phi() const {
	std::vector<TypePhi> phi(numTimepoints());
	for (size_t k = 0; k < numTimepoints(); ++k) {
		double logPhi = 0.0;
		for (size_t m = 0; m < _numEpochs; m++) { logPhi += _rho(k, m) * gamma(m); }
		phi[k] = exp(logPhi); // make logPhi to phi
	}
	return phi;
}

//------------------------------------------
// TStochasticPrior
//------------------------------------------

TStochasticPrior::TStochasticPrior(TypeParamLogPhi *LogPhi, TypeParamLogSigma *LogSigma, TypeParamGamma *Gamma,
                                   TypeParamTOC *TimesOfChange)
    : _logPhi(LogPhi), _logSigma(LogSigma), _gamma(Gamma, TimesOfChange) {}

void TStochasticPrior::initialize(const TUniqueContainer<TypeTime> &Timepoints, TBirpPrior *BirpPrior) {
	// initialize gamma and times of change
	_gamma.initialize(Timepoints, BirpPrior);

	// pre-calculate deltaTime
	_deltaTimeVec.resize(_gamma.numTimepoints(), 0);
	for (size_t k = 1; k < Timepoints.size(); k++) { _deltaTimeVec[k] = Timepoints[k] - Timepoints[k - 1]; }

	// fill rho
	_fillStochasticRho(Timepoints);

	// resize deterministic phi
	_newDeterministicPhi.resize(_gamma.numTimepoints());
	_oldDeterministicPhi.resize(_gamma.numTimepoints());

	// resize delta's
	_newSumRhoGamma.resize(_gamma.numTimepoints(), 0.0);
	_oldSumRhoGamma.resize(_gamma.numTimepoints(), 0.0);
}

void TStochasticPrior::estimateInitialGamma(const TData &Data, const TUniqueContainer<TypeTime> &Timepoints,
                                            size_t SpeciesID) {
	_gamma.estimateInitialGamma(Data, Timepoints, SpeciesID);
}

void TStochasticPrior::initializeTempVariables(const TData &Data, const TModelBase &Model, bool CalculateJP) {
	// initialize rho sum gamma
	_initializeSumRhoGamma();
	_oldSumRhoGamma = _newSumRhoGamma;

	// initialize deterministic phi
	_newDeterministicPhi = _gamma.phi();
	_oldDeterministicPhi = _newDeterministicPhi;

	// initialize JP
	if (CalculateJP) { _gamma.initializeJP(Data, Model); }

	_initializeTransitionProbabilities();
}

void TStochasticPrior::_fillStochasticRho(const TUniqueContainer<TypeTime> &Timepoints) {
	std::vector<TypeTime> epochStarts = _gamma.getEpochStarts(Timepoints);
	_stochasticRho.zeros((int)Timepoints.size(), (int)_gamma.numEpochs());

	for (size_t t = 1; t < Timepoints.size(); ++t) {
		for (size_t e = 1; e < epochStarts.size(); ++e) {
			// this is equation 6 in MS (stochastic rho)
			if (Timepoints[t] >= epochStarts[e] && Timepoints[t - 1] <= epochStarts[e - 1]) {
				_stochasticRho(t, e - 1) = epochStarts[e] - epochStarts[e - 1];
			} else if (epochStarts[e - 1] < Timepoints[t] && Timepoints[t] < epochStarts[e] &&
			           Timepoints[t - 1] <= epochStarts[e - 1]) {
				_stochasticRho(t, e - 1) = Timepoints[t] - epochStarts[e - 1];
			} else if (Timepoints[t] >= epochStarts[e] && epochStarts[e - 1] < Timepoints[t - 1] &&
			           Timepoints[t - 1] < epochStarts[e]) {
				_stochasticRho(t, e - 1) = epochStarts[e] - Timepoints[t - 1];
			} else if (epochStarts[e - 1] < Timepoints[t] && Timepoints[t] < epochStarts[e] &&
			           epochStarts[e - 1] < Timepoints[t - 1] && Timepoints[t - 1] < epochStarts[e]) {
				_stochasticRho(t, e - 1) = Timepoints[t] - Timepoints[t - 1];
			} else {
				_stochasticRho(t, e - 1) = 0.0;
			}
		}
	}
}

void TStochasticPrior::_initializeSumRhoGamma() {
	std::fill(_newSumRhoGamma.begin(), _newSumRhoGamma.end(), 0.0);
	for (size_t k = 1; k < _gamma.numTimepoints(); ++k) {
		for (size_t e = 0; e < _gamma.numEpochs(); e++) {
			_newSumRhoGamma[k] += _stochasticRho(k, e) * (double)_gamma.gamma(e);
		}
	}
}

void TStochasticPrior::_initializeTransitionProbabilities() {
	const size_t numLoc        = _logPhi->dimensions()[0];
	const size_t numTimePoints = _logPhi->dimensions()[1];

	_tryLogTransitionProbs.resize(numLoc, std::vector<double>(numTimePoints, 0.0));
	_curLogTransitionProbs.resize(numLoc, std::vector<double>(numTimePoints, 0.0));

	for (size_t j = 0; j < numLoc; ++j) {            // loop over locations
		for (size_t k = 1; k < numTimePoints; ++k) { // loop over timepoints
			_curLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
		}
	}
	_tryLogTransitionProbs = _curLogTransitionProbs;
}

coretools::probdist::TNormalDistr TStochasticPrior::_getMeanSdTransitionProbability(size_t j, size_t k) const {
	const auto deltaTime = (double)_deltaTimeVec[k];

	// calculate mean of normal distribution
	const auto sigma                = exp(_logSigma->value());
	const auto sigma2               = sigma * sigma;
	const auto sumRhoGamma          = _newSumRhoGamma[k];
	const auto log_phi_j_previous_k = _logPhi->value(_logPhi->getIndex({j, k - 1}));
	const auto mean                 = log_phi_j_previous_k + sumRhoGamma - sigma2 / 2.0 * deltaTime;

	// calculate sd of normal distribution
	const double sd = sqrt(sigma2 * deltaTime);

	coretools::probdist::TNormalDistr distr(mean, sd);
	return distr;
}

double TStochasticPrior::_calculateLogTransitionProbability(size_t j, size_t k) const {
	assert(k > 0);

	coretools::probdist::TNormalDistr distr = _getMeanSdTransitionProbability(j, k);
	const double log_phi_j_k                = _logPhi->value(_logPhi->getIndex({j, k}));
	return distr.logDensity(log_phi_j_k);
}

double TStochasticPrior::calculateLLRatioForUpdateGamma() {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) { // loop over timepoints
			_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
			LL += _tryLogTransitionProbs[j][k] - _curLogTransitionProbs[j][k];
		}
	}
	return LL;
}

double TStochasticPrior::calculateLL() const {
	// likelihood: just go through Markov chain and sum up log transition probabilities
	double LL = 0.0;
	for (size_t j = 0; j < _logPhi->dimensions()[0]; ++j) {     // loop over locations
		for (size_t k = 1; k < _logPhi->dimensions()[1]; ++k) { // loop over timepoints
			LL += _calculateLogTransitionProbability(j, k);
		}
	}
	return LL;
}

void TStochasticPrior::updateGamma(const TData &Data, const TModelBase &Model) {
	for (size_t m = 0; m < _gamma.numEpochs(); m++) { // loop over epochs and update gamma epoch by epoch
		if (_gamma.updateSpecificIndex(m)) {
			// remember current delta
			_oldSumRhoGamma = _newSumRhoGamma;
			_initializeSumRhoGamma();
			// remember current phi
			_oldDeterministicPhi       = _newDeterministicPhi;
			_newDeterministicPhi       = _gamma.phi();
			const double LLRatio       = calculateLLRatioForUpdateGamma();
			const double logPriorRatio = _gamma.getLogPriorRatio(Data, Model); // Jeffrey prior

			const double logH = LLRatio + logPriorRatio;
			if (_gamma.acceptOrReject(logH, m)) { // already swaps JP (if necessary)
				// accepted
				swapTryCur(); // -> _curLogTransitionProbs = _tryLogTransitionProbs;
			} else {
				// rejected
				_newSumRhoGamma      = _oldSumRhoGamma;
				_newDeterministicPhi = _oldDeterministicPhi;
			}
		}
	}
	// update summaries of gamma posterior after update has finished
	updateSummaryGammaPosterior();
}

bool TStochasticPrior::gammaIsUpdated() const { return _gamma.gammaIsUpdated(); }

std::vector<double> TStochasticPrior::updateGamma_forLogPhiJointUpdate(size_t m) {
	std::vector<double> Delta(_gamma.numTimepoints(), 0.0);
	if (_gamma.updateSpecificIndex(m)) {
		// remember current delta
		_oldSumRhoGamma = _newSumRhoGamma;
		_initializeSumRhoGamma();
		// remember current phi
		_oldDeterministicPhi = _newDeterministicPhi;
		_newDeterministicPhi = _gamma.phi();

		// calculate Delta
		for (size_t k = 0; k < _gamma.numTimepoints(); ++k) {
			for (size_t l = 0; l <= k; ++l) { Delta[k] += _newSumRhoGamma[l] - _oldSumRhoGamma[l]; }
		}
	}
	return Delta;
}

bool TStochasticPrior::evaluateGamma_forLogPhiJointUpdate(size_t m, double LLRatio, const TData &Data,
                                                          const TModelBase &Model) {
	const double logPriorRatio = _gamma.getLogPriorRatio(Data, Model); // Jeffrey prior
	const double logH          = LLRatio + logPriorRatio;
	if (_gamma.acceptOrReject(logH, m)) { // already swaps JP (if necessary)
		// no need to swap transition probabilities: these didn't change!
		return true;
	}
	// rejected
	_newSumRhoGamma      = _oldSumRhoGamma;
	_newDeterministicPhi = _oldDeterministicPhi;
	return false;
}

void TStochasticPrior::updateSummaryGammaPosterior() { _gamma.updateSummaryGammaPosterior(); }

double TStochasticPrior::calculateLogPriorRatio(size_t j, size_t k) {
	if (k == _gamma.numTimepoints()) { return 0.0; } // last timepoint: no data on the right

	_tryLogTransitionProbs[j][k] = _calculateLogTransitionProbability(j, k);
	return _tryLogTransitionProbs[j][k] - _curLogTransitionProbs[j][k];
}

void TStochasticPrior::swapTryCur() {
	// used when gamma updated (-> all transition probabilities change)
	std::swap(_tryLogTransitionProbs, _curLogTransitionProbs);
}

void TStochasticPrior::swapTryCur(size_t j, size_t k) {
	// used when a phi is updated
	if (k < _gamma.numTimepoints()) { // last timepoint: no data on the right
		std::swap(_tryLogTransitionProbs[j][k], _curLogTransitionProbs[j][k]);
	}
}

double TStochasticPrior::sampleFromTransitionProbability(size_t j, size_t k) const {
	coretools::probdist::TNormalDistr distr = _getMeanSdTransitionProbability(j, k);
	return distr.sample();
}

void TStochasticPrior::setSimulatedGammaToZero() { _gamma.setSimulatedGammaToZero(); }

TypePhi TStochasticPrior::deterministicPhi(size_t k) const { return _newDeterministicPhi[k]; }
std::vector<TypePhi> TStochasticPrior::calculateDeterministicPhi() const { return _gamma.phi(); }

size_t TStochasticPrior::numEpochs() const { return _gamma.numEpochs(); }

void TStochasticPrior::resetSummaryGammaPosterior() { _gamma.resetSummaryGammaPosterior(); }

void TStochasticPrior::writeSummaryGammaPosterior(const std::string &Prefix) {
	_gamma.writeSummaryGammaPosterior(Prefix);
}