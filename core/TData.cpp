//
// Created by liam on 10.03.23.
//

#include "TData.h"
#include "TEpochModel.h"
#include <regex>

//-------------------------------------------
// TZeroPoisDistr
//-------------------------------------------

class TZeroPoisDistr {
private:
	coretools::Probability _p;
	coretools::StrictlyPositive<double> _lambda;

public:
	TZeroPoisDistr(coretools::Probability p, coretools::StrictlyPositive<double> lambda);

	TZeroPoisDistr(std::string_view parameterString) { set(parameterString); }

	void set(coretools::Probability p, coretools::StrictlyPositive<double> Lambda);

	void set(std::string_view parameterString);

	static constexpr std::string_view name = "zeroPoisson";

	size_t sample(coretools::Probability p, coretools::StrictlyPositive<double> lambda) const;

	// functions using member parameters
	size_t sample() const noexcept;
};

//-------------------------------------------
// TTimepoints
//-------------------------------------------

TTimepoints::TTimepoints(const std::vector<TypeCounts> &countVec, const std::vector<TypeCovariate> &covVec,
                         const size_t &timepoint_id) {
	_counts_per_species = countVec;
	_timepoint_id       = timepoint_id;

	_covariates = covVec;

	// we (temporarily) only allow a single speces
	if (_counts_per_species.size() > 1) { UERROR("Currently only supports single species!"); }
}

TTimepoints::TTimepoints(size_t timepoint_id, size_t numSpecies, size_t numCovariates) {
	_timepoint_id = timepoint_id;
	_counts_per_species.resize(numSpecies);
	_covariates.resize(numCovariates);
}

std::vector<TypeCovariate> TTimepoints::covariates() { return _covariates; }

size_t TTimepoints::numSpecies() const { return _counts_per_species.size(); }

TypeCounts TTimepoints::counts_per_species(size_t species_id) const { return _counts_per_species[species_id]; }

TypeEffort TTimepoints::effort() const {
	return _effort;
}

TypeEffort TTimepoints::effort_old() { return _effort_old; }

size_t TTimepoints::timepoint_id() const { return _timepoint_id; }

void TTimepoints::setTimepoint(size_t timepointToAdd) { _timepoint_id = timepointToAdd; }

void TTimepoints::setCovariate(size_t c, TypeCovariate Covariate) { _covariates[c] = Covariate; }

void TTimepoints::initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	if (_covariates.size() == 1) {
		_effort = _covariates[0];
	} else {
		// calculate s_{ijk} = sum_{c=1}^C x_{ijkc} * beta_{ic}
		_effort = 0.0;
		for (size_t c = 0; c < _covariates.size(); ++c) { _effort += _covariates[c] * Beta->value(c); }
	}
	_effort_old = _effort;
}

void TTimepoints::updateEffort(TypeBeta NewBeta, TypeBeta OldBeta, size_t c) {
	// remember current effort
	_effort_old = _effort;
	// update effort: no need to re-calculate sum over all beta, just update with difference
	// s_{ijk}' = s_{ijk} + x_{ijkc} * (beta_{ic}' - beta_{ic})
	_effort += _covariates[c] * (NewBeta - OldBeta);
}

void TTimepoints::resetEffort() {
	// this function is called in case beta update is rejected
	_effort = _effort_old;
}

void TTimepoints::write(coretools::TOutputFile &File, std::string location,
                        const TUniqueContainer<size_t> &Timepoints) const {
	File << location << Timepoints[_timepoint_id] << _counts_per_species << _covariates << coretools::endl;
}

void TTimepoints::_simulatePoisson(double lambda) {
	lambda *= _effort;
	if (lambda == 0.0) {
		_counts_per_species[0] = 0;
	} else {
		_counts_per_species[0] = coretools::instances::randomGenerator().getPoissonRandom((double)lambda);
	}
}

void TTimepoints::simulatePoisson(size_t N0, double phi, double detectionProbs) {
	coretools::Positive<double> lambda = (double)N0 * phi * detectionProbs;
	_simulatePoisson(lambda);
}

void TTimepoints::simulatePoisson(double phi) { _simulatePoisson(phi); }

void TTimepoints::simulateNB(double phi, double a, double mu_ij) {
	// first draw lambda from a Gamma distr and then from there draw counts from Pois
	double alpha           = a * mu_ij * phi * _effort;
	double beta            = a;
	double lambda          = coretools::instances::randomGenerator().getGammaRand(alpha, beta);
	_counts_per_species[0] = coretools::instances::randomGenerator().getPoissonRandom(lambda);
}

//-------------------------------------------
// TLocations
//-------------------------------------------

size_t TLocations::size() const { return _timepoints.size(); }

TypeNu TLocations::nu(size_t species_id) const { return _nu[species_id]; }

size_t TLocations::location_id() const { return _location_id; }

TLocations::TLocations(const TTimepoints &pt, size_t location_id) {
	_location_id = location_id;
	add(pt);
}

TLocations::TLocations(size_t location_id, size_t TimepointSize, size_t numSpecies, size_t numCovariates) {
	_location_id = location_id;
	for (size_t t = 0; t < TimepointSize; t++) { _timepoints.emplace_back(t, numSpecies, numCovariates); }
}

const TTimepoints &TLocations::operator[](const size_t index) const { return _timepoints[index]; }

void TLocations::add(const TTimepoints &pt) {
	for (auto &t : _timepoints) {
		if (t.timepoint_id() == pt.timepoint_id()) {
			UERROR("Make sure your ipnut file does not contain multiple rows with identical location and timepoints");
		}
	}
	// when we are here: timpoint is unique -> push it back
	_timepoints.push_back(pt);
	// now add nu, first check if nu exists, otherwise create it
	if (_nu.empty()) { _nu.resize(pt.numSpecies(), 0); }
	for (size_t i = 0; i < pt.numSpecies(); i++) { _nu[i] += (uint32_t)pt.counts_per_species(i); }
}

void TLocations::checkEffortBiggerThanZero() {}

void TLocations::sorttimes(const std::vector<size_t> &index) {
	std::vector<size_t> positionOfTimepointIds;
	for (auto &timepoint : _timepoints) {
		auto it    = std::find(index.begin(), index.end(), timepoint.timepoint_id());
		size_t pos = it - index.begin();
		positionOfTimepointIds.push_back(pos);
	}
	// sort the newly found positions of timepoints using ranksort()
	std::vector<size_t> sortingIndex = coretools::rankSort(positionOfTimepointIds);
	// sort _timepoints
	_timepoints                      = coretools::sortContainerByRank(_timepoints, sortingIndex);
	// now sort positionOfTimepointIds
	positionOfTimepointIds           = coretools::sortContainerByRank(positionOfTimepointIds, sortingIndex);
	// now set timepoint_id
	for (size_t i = 0; i < _timepoints.size(); i++) { _timepoints[i].setTimepoint(positionOfTimepointIds[i]); }
}

void TLocations::write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
                       const TUniqueContainer<size_t> &Timepoints) const {
	std::string location = Locations[_location_id];
	for (const auto &timepoint : _timepoints) { timepoint.write(File, location, Timepoints); }
}

std::vector<TTimepoints>::const_iterator TLocations::cbegin() const { return _timepoints.cbegin(); }

std::vector<TTimepoints>::iterator TLocations::begin() { return _timepoints.begin(); }

std::vector<TTimepoints>::const_iterator TLocations::cend() const { return _timepoints.cend(); }

std::vector<TTimepoints>::iterator TLocations::end() { return _timepoints.end(); }

void TLocations::initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	for (auto &timepoint : _timepoints) { timepoint.initializeEffort(Beta); }
}

void TLocations::updateEffort(TypeBeta NewBeta, TypeBeta OldBeta, size_t c) {
	for (auto &timepoint : _timepoints) { timepoint.updateEffort(NewBeta, OldBeta, c); }
}

void TLocations::resetEffort() {
	for (auto &timepoint : _timepoints) { timepoint.resetEffort(); }
}

void TLocations::simulatePoisson(size_t N0, const std::vector<double> &phi, double detectionProbs) {
	for (size_t t = 0; t < _timepoints.size(); t++) { _timepoints[t].simulatePoisson(N0, phi[t], detectionProbs); }
}

void TLocations::simulatePoisson(const std::vector<double> &phi) {
	for (size_t t = 0; t < _timepoints.size(); t++) { _timepoints[t].simulatePoisson(phi[t]); }
}

void TLocations::simulateNB(const std::vector<double> &phi, double a, double mu_ij) {
	for (size_t t = 0; t < _timepoints.size(); t++) { _timepoints[t].simulateNB(phi[t], a, mu_ij); }
}

//-------------------------------------------
// TMethods
//-------------------------------------------

TMethods::TMethods(const std::string_view &methodName, size_t index,
                   const std::vector<std::size_t> &speciesIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateIDsinUniqueContainer)
    : _methodName(methodName), _method_id(index), _covariateIDsinUniqueContainer(covariateIDsinUniqueContainer),
      _speciesIDsinUniqueContainer(speciesIDsinUniqueContainer) { // infer or simulate from file
}

TMethods::TMethods(const std::string &methodName, size_t locationSize, size_t TimepointSize, size_t index,
                   const std::vector<size_t> &speciesIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateIDsinUniqueContainer)
    : _methodName(methodName), _method_id(index), _covariateIDsinUniqueContainer(covariateIDsinUniqueContainer),
      _speciesIDsinUniqueContainer(speciesIDsinUniqueContainer) { // simulate

	_locations.resize(locationSize);
	for (size_t i = 0; i < locationSize; i++) {
		_locations.emplace_back(i, TimepointSize, _speciesIDsinUniqueContainer.size(),
		                        covariateIDsinUniqueContainer.size());
	}
}

void TMethods::addBeta(std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) { _beta = Beta; }

void TMethods::add(const TTimepoints &pt, size_t location_id) {
	bool dataPointExists = false;
	for (auto &loc : _locations) {
		if (loc.location_id() == location_id) {
			dataPointExists = true;
			loc.add(pt);
			break;
		}
	}
	if (!dataPointExists) { _locations.emplace_back(pt, location_id); }
}

void TMethods::fillIDVectors(size_t totalSizeOfUniqueContainer) {
	_hasDataForSpeciesID.resize(totalSizeOfUniqueContainer, false);
	_indexInCounts.resize(totalSizeOfUniqueContainer, 0);
	for (size_t c = 0; c < _speciesIDsinUniqueContainer.size(); c++) {
		_hasDataForSpeciesID[_speciesIDsinUniqueContainer[c]] = true;
		_indexInCounts[_speciesIDsinUniqueContainer[c]]       = c;
	}
}

void TMethods::checkEffortBiggerThanZero() { // todo: check if used somewhere, else use it :-)
	for (auto &loc : _locations) { loc.checkEffortBiggerThanZero(); }
}

void TMethods::sorttimes(const std::vector<size_t> &sortingIndex) {
	for (auto &loc : _locations) { loc.sorttimes(sortingIndex); }
}

void TMethods::initializeEffort() {
	// initialize effort
	for (auto &location : _locations) { location.initializeEffort(_beta); }
}

size_t TMethods::size() const { return _locations.size(); }
const std::string &TMethods::name() const { return _methodName; }

size_t TMethods::method_id() const { return _method_id; }
size_t TMethods::numCovariates() const { return _covariateIDsinUniqueContainer.size(); }

bool TMethods::hasDataForSpeciesID(size_t speciesIDinUniqueContainer) const {
	return _hasDataForSpeciesID[speciesIDinUniqueContainer];
}

size_t TMethods::getIndexInCounts(size_t speciesIDinUniqueContainer) const {
	return _indexInCounts[speciesIDinUniqueContainer];
}

const std::vector<size_t> &TMethods::speciesIDsinUniqueContainer() const { return _speciesIDsinUniqueContainer; }
const std::vector<size_t> &TMethods::covariateIDsinUniqueContainer() const { return _covariateIDsinUniqueContainer; }
size_t TMethods::covariateIDinUniqueContainer(size_t c) const { return _covariateIDsinUniqueContainer[c]; }

std::vector<TLocations>::const_iterator TMethods::cbegin() const { return _locations.cbegin(); }

std::vector<TLocations>::iterator TMethods::begin() { return _locations.begin(); }

std::vector<TLocations>::const_iterator TMethods::cend() const { return _locations.cend(); }

std::vector<TLocations>::iterator TMethods::end() { return _locations.end(); }

std::vector<std::string> TMethods::_splitVecWithDistributions(const std::string &Input) {
	std::vector<std::string> output;

	std::regex regex(",(?![^()]*\\))"); // matches comma not inside brackets
	std::sregex_token_iterator it(Input.begin(), Input.end(), regex, -1);
	std::sregex_token_iterator end;
	for (; it != end; ++it) { output.push_back(*it); }

	return output;
}

void TMethods::_simulateOneCovariate(size_t c, const std::string &Cov) {
	if (coretools::str::stringIsProbablyANumber(Cov)) {
		// case 1: is a number: use it for all covariates
		double value = coretools::str::fromString<double>(Cov);
		for (auto &_location : _locations) {
			for (auto &t : _location) { t.setCovariate(c, value); }
		}
	} else {
		// case 2: parse distributions. "Poisson(1), uniform(a,b), zero-poisson(p, lambda) -> p: p(not 0),
		coretools::str::TSplitter sp(Cov, '(');
		std::string_view distr = sp.front();          // takes first element
		sp.popFront();                                // takes 'Poisson' away
		std::string params = (std::string)sp.front(); // now parameters are in front
		params             = coretools::str::extractBefore(params, ')');
		if (distr == "Poisson") {
			_simulateCovariateFromDistribution<coretools::probdist::TPoissonDistr>(c, params);
		} else if (distr == "uniform") {
			_simulateCovariateFromDistribution<coretools::probdist::TUniformDistr>(c, params);
		} else if (distr == "zeroPoisson") {
			_simulateCovariateFromDistribution<TZeroPoisDistr>(c, params);
		} else {
			UERROR("Unknown distribution '", distr,
			       "' (argument 'covariates'). Supported distributions are: Poisson, uniform, zeroPoisson.");
		}
	}
}

void TMethods::_simulateCovariate() {
	using namespace coretools;
	std::string effortString = instances::parameters().getParameterWithDefault("covariates", "Poisson(1000)");
	if (str::stringIsProbablyANumber(effortString)) {
		// case 1: is a number: use it for all covariates
		for (size_t c = 0; c < numCovariates(); ++c) { _simulateOneCovariate(c, effortString); }
	} else if (str::stringContains(effortString, ',')) {
		// cases 2 and 4: comma-separated string of numbers (case 2) or distributions (case 4) or both
		std::vector<std::string> vec = _splitVecWithDistributions(effortString);
		if (vec.size() != numCovariates()) {
			UERROR("The size of covariates provided with argument 'covariates' (", vec.size(),
			       ") does not match the number of covariates (", numCovariates(), ").");
		}
		for (size_t c = 0; c < vec.size(); ++c) { _simulateOneCovariate(c, vec[c]); }
	} else {
		// case 3: use the same distribution for all covariates
		for (size_t c = 0; c < numCovariates(); ++c) { _simulateOneCovariate(c, effortString); }
	}
	// now that we have simulated all covariates: re-calculate effort
	initializeEffort();
}

// function to simulate
std::pair<size_t, size_t> TMethods::_getIndexInDetectionProbabilitiesToSimulate(
    const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
    const TUniqueContainer<std::string> &locations, size_t i) {
	std::string locName  = locations[_locations[i].location_id()];
	std::string methName = _methodName;
	size_t locIndex      = detectionProbs.getDimensionName(0)->getIndex(locName);
	size_t methIndex     = detectionProbs.getDimensionName(1)->getIndex(methName);
	return {locIndex, methIndex};
}

void TMethods::simulatePoisson(size_t N0, const std::vector<double> &phi,
                               const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
                               const TUniqueContainer<std::string> &locations) {
	for (size_t i = 0; i < size(); i++) {
		auto [locIndexInDetection, methIndexInDetection] =
		    _getIndexInDetectionProbabilitiesToSimulate(detectionProbs, locations, i);
		_locations[i].simulatePoisson(N0, phi, detectionProbs(locIndexInDetection, methIndexInDetection));
	}
}

void TMethods::simulatePoisson(const std::vector<double> &phi) {
	_simulateCovariate();
	for (size_t i = 0; i < size(); i++) { _locations[i].simulatePoisson(phi); }
}

std::vector<double> TMethods::simulateNB(size_t N0, const std::vector<double> &phi,
                                         const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
                                         const TUniqueContainer<std::string> &locations, double A) {
	std::vector<double> mus(_locations.size());
	for (size_t i = 0; i < size(); i++) {
		auto [locIndexInDetection, methIndexInDetection] =
		    _getIndexInDetectionProbabilitiesToSimulate(detectionProbs, locations, i);
		mus[i] = (double)N0 * detectionProbs(locIndexInDetection, methIndexInDetection);
		_locations[i].simulateNB(phi, A, mus[i]);
	}
	return mus;
}

void TMethods::simulateNB(const std::vector<double> &phi, double a, const std::vector<double> &mu,
                          bool SimulateEffort) {
	if (SimulateEffort) { _simulateCovariate(); }
	for (size_t i = 0; i < size(); i++) { _locations[i].simulateNB(phi, a, mu[i]); }
}

void TMethods::write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
                     const TUniqueContainer<size_t> &Timepoints) const {
	for (const auto &location : _locations) { location.write(File, Locations, Timepoints); }
}

void TMethods::initializeInferred(const TUniqueContainer<std::string> &CovariateNames) {
	// initialize beta: first get appropriate covariate name
	auto covNames = std::make_shared<coretools::TNamesStrings>();
	for (size_t c = 0; c < numCovariates(); ++c) {
		const size_t id = _covariateIDsinUniqueContainer[c];
		covNames->addName({CovariateNames[id]});
	}
	_beta->initStorage({numCovariates()}, {covNames});
}

void TMethods::estimateInitialBeta() {
	// TODO: How to initialize beta?
	// set initial value of all beta's to one
	for (size_t c = 0; c < _beta->totalSize(); ++c) { _beta->setInitial(c, 1.0); }
	// set effort based on initial values of beta
	initializeEffort();
	if (numCovariates() == 1) { // if only a single covariate: do not update beta
		_beta->setIsUpdated(false);
	}
}

void TMethods::_updateBeta(size_t c, const std::vector<std::unique_ptr<TEpochModel>> &EpochModels) {
	// propose new value
	if (c == 1) { return; }
	if (!_beta->updateSpecificIndex(c)) { return; } // beta is not updated -> stop here
	// update effort for all locations and timepoints of current method
	for (auto &location : _locations) { location.updateEffort(_beta->value(c), _beta->oldValue(c), c); }
	// calculate Hastings ratio
	double logRatio = _beta->getLogPriorRatio(c);
	std::vector<std::vector<double>> LLs(EpochModels.size());
	for (size_t s = 0; s < EpochModels.size(); ++s) {
		LLs[s] = EpochModels[s]->calculateLLForBetaUpdate(*this);
		logRatio += coretools::containerSum(LLs[s]) - coretools::containerSum(EpochModels[s]->oldLLs());
	}
	// decide if accept or reject
	if (_beta->acceptOrReject(logRatio, c)) {
		// accepted: set likelihood in all epoch models
		for (size_t s = 0; s < EpochModels.size(); ++s) { EpochModels[s]->setLLs(LLs[s]); }
	} else {
		// rejected: re-set effort to old value before update
		for (auto &location : _locations) { location.resetEffort(); }
	}
}

void TMethods::updateBeta(const std::vector<std::unique_ptr<TEpochModel>> &EpochModels) {
	for (size_t c = 0; c < numCovariates(); ++c) { _updateBeta(c, EpochModels); }
}

//-------------------------------------------
// TZeroPoisDistr
//-------------------------------------------

TZeroPoisDistr::TZeroPoisDistr(coretools::Probability p, coretools::StrictlyPositive<double> lambda)
    : _p(p), _lambda(lambda) {}

void TZeroPoisDistr::set(coretools::Probability p, coretools::StrictlyPositive<double> Lambda) {
	_lambda = Lambda;
	_p      = p;
}

void TZeroPoisDistr::set(std::string_view parameterString) {
	coretools::str::convertString(
	    parameterString, std::string{"Use "}.append(name) + "(p, lambda) with 0 <= p <= 1, lambda>0.", _p, _lambda);
}

size_t TZeroPoisDistr::sample(coretools::Probability p, coretools::StrictlyPositive<double> lambda) const {
	double tmp = coretools::instances::randomGenerator().getRand(0, 1);
	if (tmp > p) { // p is probability of taking Poisson
		return 0;
	} else {
		return coretools::instances::randomGenerator().getPoissonRandom(lambda);
	}
}

size_t TZeroPoisDistr::sample() const noexcept { return TZeroPoisDistr::sample(_p, _lambda); }