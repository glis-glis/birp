//
// Created by liam on 10.03.23.
//

#include "TData.h"
#include "TModels.h"
#include "stattools/ParametersObservations/TIndexSampler.h"
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

TTimepoints::TTimepoints(const std::vector<TypeCounts> &Counts, const std::vector<TypeCovariateEffort> &CovsEffort,
                         const std::vector<TypeCovariateDetection> &CovsDetection, const size_t &Timepoint_id)
    : _counts_per_species(Counts), _covariatesEffort(CovsEffort), _covariatesDetection(CovsDetection),
      _timepoint_id(Timepoint_id) {

	// we (temporarily) only allow a single speces
	if (_counts_per_species.size() > 1) { UERROR("Currently only supports single species!"); }
}

TTimepoints::TTimepoints(size_t timepoint_id, size_t numSpecies, size_t numCovariatesEffort,
                         size_t numCovariatesDetection) {
	_timepoint_id = timepoint_id;
	_counts_per_species.resize(numSpecies);
	_covariatesEffort.resize(numCovariatesEffort);
	_covariatesDetection.resize(numCovariatesDetection);
}

const std::vector<TypeCovariateEffort> &TTimepoints::covariatesEffort() const { return _covariatesEffort; }

const std::vector<TypeCovariateDetection> &TTimepoints::covariatesDetection() const { return _covariatesDetection; }

size_t TTimepoints::numSpecies() const { return _counts_per_species.size(); }

TypeCounts TTimepoints::counts_per_species(size_t species_id) const { return _counts_per_species[species_id]; }

TypeEffort TTimepoints::effort() const { return _effort; }

TypeEffort TTimepoints::effort_old() const { return _effort_old; }

size_t TTimepoints::timepoint_id() const { return _timepoint_id; }

void TTimepoints::setTimepoint(size_t timepointToAdd) { _timepoint_id = timepointToAdd; }

void TTimepoints::setCovariateEffort(size_t c, TypeCovariateEffort Covariate) { _covariatesEffort[c] = Covariate; }
void TTimepoints::setCovariateDetection(size_t c, TypeCovariateDetection Covariate) {
	_covariatesDetection[c] = Covariate;
}

bool TTimepoints::inferAlpha() const {
	// do not infer if single covariate for effort
	// TODO: confirm with Daniel
	return _covariatesEffort.size() > 1;
}

bool TTimepoints::inferBeta() const {
	// do not infer if single covariate for effort and no covariates for detection
	// TODO: confirm with Daniel
	return _covariatesDetection.size() > 1;
}

void TTimepoints::_standardizeCovariatesDetection(double MeanDetectionCovariates, double SdDetectionCovariates) {
	if (_covariatesDetection.size() < 2) { return; } // zero or one detection covariates: no need to do anything

	// standardize detection probability covariates with standard score (mean zero, unit variance)
	for (auto &cov : _covariatesDetection) { cov = (cov - MeanDetectionCovariates) / SdDetectionCovariates; }
}

void TTimepoints::_standardizeCovariatesEffort() {
	assert(!_covariatesEffort.empty()); // should always have at least one effort covariate

	// standardize effort covariates such that they sum to one
	coretools::normalize(_covariatesEffort);
}

void TTimepoints::standardizeCovariates(double MeanDetectionCovariates, double SdDetectionCovariates) {
	_standardizeCovariatesDetection(MeanDetectionCovariates, SdDetectionCovariates);
	_standardizeCovariatesEffort();
}

double TTimepoints::_calculateEffortTerm(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha) {
	assert(!_covariatesEffort.empty());
	if (!inferAlpha()) { return _covariatesEffort[0]; }

	// calculate sum_{c=1}^C alpha_{ic} x_{ijkc}
	double effortTerm = 0.0;
	for (size_t c = 0; c < _covariatesEffort.size(); ++c) { effortTerm += Alpha->value(c) * _covariatesEffort[c]; }
	return effortTerm;
}

double TTimepoints::_calculateDetectionTerm(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                                            const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	if (_covariatesDetection.empty()) { return 1.0; }
	if (!inferBeta()) { return coretools::logistic(_covariatesDetection[0]); }

	// calculate logistic(beta_{i0} + sum_{d=1}^D * beta_{id} y_{ijkd})
	double detectionTerm = Beta0->value();
	for (size_t d = 0; d < _covariatesDetection.size(); ++d) {
		detectionTerm += Beta->value(d) * _covariatesDetection[d];
	}
	return coretools::logistic(detectionTerm);
}

void TTimepoints::initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
                                   const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                                   const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	// update alpha and beta
	// calculate s_{ijk} = sum_{c=1}^C alpha_{ic} x_{ijkc} * logistic(beta_{i0} + sum_{d=1}^D * beta_{id} y_{ijkd})
	_effortTerm    = _calculateEffortTerm(Alpha);
	_detectionTerm = _calculateDetectionTerm(Beta0, Beta);
	_effort        = _effortTerm * _detectionTerm;

	// also set old
	_effortTerm_old    = _effortTerm;
	_detectionTerm_old = _detectionTerm;
	_effort_old        = _effort;
}

void TTimepoints::updateAlphaEffort(TypeAlpha NewAlpha_c1, TypeAlpha OldAlpha_c1, TypeAlpha NewAlpha_c2,
                                    TypeAlpha OldAlpha_c2, size_t c1, size_t c2) {
	// remember current effort
	_effort_old     = _effort;
	_effortTerm_old = _effortTerm;

	// update effort: no need to re-calculate sum over all alpha, just update with difference for both alpha
	// s_{ijk}' = s_{ijk} + detection_term * x_{ijkc} * (alpha_{ic}' - alpha_{ic})
	_effortTerm += _covariatesEffort[c1] * ((double)NewAlpha_c1 - (double)OldAlpha_c1) +
	               _covariatesEffort[c2] * ((double)NewAlpha_c2 - (double)OldAlpha_c2);
	_effort = _effortTerm * _detectionTerm;
}

void TTimepoints::updateBetaEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                                   const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	// remember current effort
	_effort_old        = _effort;
	_detectionTerm_old = _detectionTerm;

	_detectionTerm = _calculateDetectionTerm(Beta0, Beta);
	_effort        = _effortTerm * _detectionTerm;
}

void TTimepoints::resetAlphaEffort() {
	// this function is called in case alpha update is rejected
	_effortTerm = _effortTerm_old;
	_effort     = _effort_old;
}

void TTimepoints::resetBetaEffort() {
	// this function is called in case beta update is rejected
	_detectionTerm = _detectionTerm_old;
	_effort        = _effort_old;
}

void TTimepoints::write(coretools::TOutputFile &File, std::string location,
                        const TUniqueContainer<size_t> &Timepoints) const {
	File << location << Timepoints[_timepoint_id] << _counts_per_species;

	// write effort covariate(s)
	for (const auto &it : _covariatesEffort) { File << it; }

	// write detection covariate(s)
	if (_covariatesDetection.size() == 1) {
		File << _covariatesDetection[0];
	} else if (_covariatesDetection.size() > 1) {
		// also contains an intercept, but this one should not be written!
		for (auto it = _covariatesDetection.begin() + 1; it != _covariatesDetection.end(); ++it) { File << *it; }
	}

	File << coretools::endl;
}

void TTimepoints::simulatePoisson(size_t N0, double phi, double detectionProbs) {
	coretools::Positive<double> lambda = (double)N0 * phi * detectionProbs;

	lambda *= _effort;
	if (lambda == 0.0) {
		_counts_per_species[0] = 0;
	} else {
		_counts_per_species[0] = coretools::instances::randomGenerator().getPoissonRandom((double)lambda);
	}
}

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

TLocations::TLocations(const TTimepoints &pt, size_t location_id) {
	_location_id = location_id;
	add(pt);
}

TLocations::TLocations(size_t location_id, size_t TimepointSize, size_t numSpecies, size_t numCovariatesEffort,
                       size_t numCovariatesDetection) {
	_location_id = location_id;
	for (size_t t = 0; t < TimepointSize; t++) {
		_timepoints.emplace_back(t, numSpecies, numCovariatesEffort, numCovariatesDetection);
	}
}

size_t TLocations::size() const { return _timepoints.size(); }

TypeNu TLocations::nu(size_t species_id) const { return _nu[species_id]; }

size_t TLocations::location_id() const { return _location_id; }

const TTimepoints &TLocations::operator[](const size_t index) const { return _timepoints[index]; }

void TLocations::add(const TTimepoints &pt) {
	for (auto &t : _timepoints) {
		if (t.timepoint_id() == pt.timepoint_id()) {
			UERROR("Make sure your input file does not contain multiple rows with identical location and timepoints");
		}
	}
	// when we are here: timpoint is unique -> push it back
	_timepoints.push_back(pt);
	// now add nu, first check if nu exists, otherwise create it
	if (_nu.empty()) { _nu.resize(pt.numSpecies(), 0); }
	for (size_t i = 0; i < pt.numSpecies(); i++) { _nu[i] += (uint64_t)pt.counts_per_species(i); }
}

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
	const std::string &location = Locations[_location_id];
	for (const auto &timepoint : _timepoints) { timepoint.write(File, location, Timepoints); }
}

std::vector<TTimepoints>::const_iterator TLocations::cbegin() const { return _timepoints.cbegin(); }

std::vector<TTimepoints>::iterator TLocations::begin() { return _timepoints.begin(); }

std::vector<TTimepoints>::const_iterator TLocations::cend() const { return _timepoints.cend(); }

std::vector<TTimepoints>::iterator TLocations::end() { return _timepoints.end(); }

void TLocations::initializeEffort(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
                                  const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                                  const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	for (auto &timepoint : _timepoints) { timepoint.initializeEffort(Alpha, Beta0, Beta); }
}

void TLocations::updateAlphaEffort(TypeAlpha NewAlpha_c1, TypeAlpha OldAlpha_c1, TypeAlpha NewAlpha_c2,
                                   TypeAlpha OldAlpha_c2, size_t c1, size_t c2) {
	for (auto &timepoint : _timepoints) {
		timepoint.updateAlphaEffort(NewAlpha_c1, OldAlpha_c1, NewAlpha_c2, OldAlpha_c2, c1, c2);
	}
}

void TLocations::updateBetaEffort(const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                                  const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	for (auto &timepoint : _timepoints) { timepoint.updateBetaEffort(Beta0, Beta); }
}

void TLocations::resetAlphaEffort() {
	for (auto &timepoint : _timepoints) { timepoint.resetAlphaEffort(); }
}

void TLocations::resetBetaEffort() {
	for (auto &timepoint : _timepoints) { timepoint.resetBetaEffort(); }
}

void TLocations::simulatePoisson(size_t N0, const TModelBase &Model, double detectionProbs) {
	for (auto &timepoint : _timepoints) {
		timepoint.simulatePoisson(N0, Model.phi(location_id(), timepoint.timepoint_id()), detectionProbs);
	}
}

void TLocations::simulateNB(const TModelBase &Model, double a, double mu_ij) {
	for (auto &timepoint : _timepoints) {
		timepoint.simulateNB(Model.phi(location_id(), timepoint.timepoint_id()), a, mu_ij);
	}
}

//-------------------------------------------
// TMethods
//-------------------------------------------

TMethods::TMethods(const std::string_view &methodName, size_t index,
                   const std::vector<std::size_t> &speciesIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateEffortIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateDetectionIDsinUniqueContainer)
    : _methodName(methodName), _method_id(index),
      _covariateEffortIDsinUniqueContainer(covariateEffortIDsinUniqueContainer),
      _covariateDetectionIDsinUniqueContainer(covariateDetectionIDsinUniqueContainer),
      _speciesIDsinUniqueContainer(speciesIDsinUniqueContainer) { // infer or simulate from file
}

TMethods::TMethods(const std::string &methodName, size_t locationSize, size_t TimepointSize, size_t index,
                   const std::vector<size_t> &speciesIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateEffortIDsinUniqueContainer,
                   const std::vector<std::size_t> &covariateDetectionIDsinUniqueContainer)
    : _methodName(methodName), _method_id(index),
      _covariateEffortIDsinUniqueContainer(covariateEffortIDsinUniqueContainer),
      _covariateDetectionIDsinUniqueContainer(covariateDetectionIDsinUniqueContainer),
      _speciesIDsinUniqueContainer(speciesIDsinUniqueContainer) { // simulate

	_locations.reserve(locationSize);
	for (size_t i = 0; i < locationSize; i++) {
		_locations.emplace_back(i, TimepointSize, _speciesIDsinUniqueContainer.size(),
		                        covariateEffortIDsinUniqueContainer.size(),
		                        covariateDetectionIDsinUniqueContainer.size());
	}
}

void TMethods::addAlphaBeta(const std::shared_ptr<stattools::TParameterTyped<TypeAlpha, 1>> &Alpha,
                            const std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta0,
                            std::shared_ptr<stattools::TParameterTyped<TypeBeta, 1>> &Beta) {
	_alpha = Alpha;
	_beta0 = Beta0;
	_beta  = Beta;
}

auto getCovNames(const TUniqueContainer<std::string> &CovariateNames, const std::vector<size_t> &ID) {
	auto covNames = std::make_shared<coretools::TNamesStrings>();
	for (size_t c = 0; c < CovariateNames.size(); ++c) {
		const size_t id = ID[c];
		covNames->addName({CovariateNames[id]});
	}
	return covNames;
}

void TMethods::initialize(const TUniqueContainer<std::string> &CovariateEffortNames,
                          const TUniqueContainer<std::string> &CovariateDetectionNames) {
	// initialize alpha: get names from effort covariates
	auto covNamesAlpha = getCovNames(CovariateEffortNames, _covariateEffortIDsinUniqueContainer);
	_alpha->initStorage({numCovariatesEffort()}, {covNamesAlpha});

	// initialize beta0: single element
	_beta0->initStorage();

	// initialize beta: get names from detection covariates
	auto covNamesBeta = getCovNames(CovariateDetectionNames, _covariateDetectionIDsinUniqueContainer);
	_beta->initStorage({numCovariatesEffort()}, {covNamesBeta});
}

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

void TMethods::standardizeCovariates() {
	// calculate mean and variance across all detection covariates
	// but check per detection covariate for unique values (not allowed)
	coretools::TMeanVar<double> meanVar;
	// TODO: check for unique covariates per column or across all?
	std::vector<std::set<double>> uniqueValues_perCovariate(_covariateDetectionIDsinUniqueContainer.size());
	for (auto &loc : _locations) {
		for (auto &tp : loc) {
			const auto &covariates = tp.covariatesDetection();
			for (size_t d = 0; d < covariates.size(); ++d) {
				meanVar.add(covariates[d]);
				uniqueValues_perCovariate[d].insert(covariates[d]);
			}
		}
	}

	// check if all covariates of all detection are unique
	for (const auto &d : uniqueValues_perCovariate) {
		if (d.size() == 1) {
			UERROR("All detection covariates are equal (", *d.begin(),
			       "). Please exclude this covariate or provide non-equal values");
		}
	}

	// now standardize
	for (auto &loc : _locations) {
		for (auto &timepoint : loc) { timepoint.standardizeCovariates(meanVar.mean(), meanVar.sd()); }
	}
}

void TMethods::fillIDVectors(size_t totalSizeOfUniqueContainer) {
	_hasDataForSpeciesID.resize(totalSizeOfUniqueContainer, false);
	_indexInCounts.resize(totalSizeOfUniqueContainer, 0);
	for (size_t c = 0; c < _speciesIDsinUniqueContainer.size(); c++) {
		_hasDataForSpeciesID[_speciesIDsinUniqueContainer[c]] = true;
		_indexInCounts[_speciesIDsinUniqueContainer[c]]       = c;
	}
}

void TMethods::sorttimes(const std::vector<size_t> &sortingIndex) {
	for (auto &loc : _locations) { loc.sorttimes(sortingIndex); }
}

void TMethods::initializeEffort() {
	_alpha->normalizeToSum1AfterUpdate();

	// initialize effort
	for (auto &location : _locations) { location.initializeEffort(_alpha, _beta0, _beta); }
}

size_t TMethods::size() const { return _locations.size(); }
const std::string &TMethods::name() const { return _methodName; }

size_t TMethods::method_id() const { return _method_id; }
size_t TMethods::numCovariatesEffort() const { return _covariateEffortIDsinUniqueContainer.size(); }
size_t TMethods::numCovariatesDetection() const { return _covariateDetectionIDsinUniqueContainer.size(); }
bool TMethods::inferAlpha() const { return _locations.cbegin()->cbegin()->inferAlpha(); }
bool TMethods::inferBeta() const { return _locations.cbegin()->cbegin()->inferBeta(); }

bool TMethods::hasDataForSpeciesID(size_t speciesIDinUniqueContainer) const {
	return _hasDataForSpeciesID[speciesIDinUniqueContainer];
}

size_t TMethods::getIndexInCounts(size_t speciesIDinUniqueContainer) const {
	return _indexInCounts[speciesIDinUniqueContainer];
}

const std::vector<size_t> &TMethods::speciesIDsinUniqueContainer() const { return _speciesIDsinUniqueContainer; }
const std::vector<size_t> &TMethods::covariateEffortIDsinUniqueContainer() const {
	return _covariateEffortIDsinUniqueContainer;
}
size_t TMethods::covariateEffortIDinUniqueContainer(size_t c) const { return _covariateEffortIDsinUniqueContainer[c]; }
const std::vector<size_t> &TMethods::covariateDetectionIDsinUniqueContainer() const {
	return _covariateDetectionIDsinUniqueContainer;
}
size_t TMethods::covariateDetectionIDinUniqueContainer(size_t c) const {
	return _covariateDetectionIDsinUniqueContainer[c];
}

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

std::pair<std::string, std::string> TMethods::_getParamsAndDistr(const std::string &Cov) {
	coretools::str::TSplitter sp(Cov, '(');
	std::string_view distr = sp.front();          // takes first element
	sp.popFront();                                // takes 'Poisson' away
	std::string params = (std::string)sp.front(); // now parameters are in front
	params             = coretools::str::extractBefore(params, ')');
	return {(std::string)distr, params};
}

void TMethods::_simulateOneCovariateEffort(size_t c, const std::string &Cov) {
	if (coretools::str::stringIsProbablyANumber(Cov)) {
		// case 1: is a number: use it for all covariates
		double value = coretools::str::fromString<double>(Cov);
		for (auto &_location : _locations) {
			for (auto &t : _location) { t.setCovariateEffort(c, value); }
		}
	} else {
		// case 2: parse distributions. "Poisson(1), uniform(a,b), zero-poisson(p, lambda) -> p: p(not 0),
		const auto [distr, params] = _getParamsAndDistr(Cov);
		if (distr == "Poisson") {
			_simulateCovariateFromDistribution<coretools::probdist::TPoissonDistr, true>(c, params);
		} else if (distr == "uniform") {
			_simulateCovariateFromDistribution<coretools::probdist::TUniformDistr, true>(c, params);
		} else if (distr == "zeroPoisson") {
			_simulateCovariateFromDistribution<TZeroPoisDistr, true>(c, params);
		} else {
			UERROR("Unknown distribution '", distr,
			       "' (argument 'covariatesEffort'). Supported distributions are: Poisson, uniform, zeroPoisson.");
		}
	}
}

void TMethods::_simulateCovariatesEffort() {
	using namespace coretools;
	std::string effortString = instances::parameters().getParameterWithDefault("covariatesEffort", "Poisson(1000)");
	if (str::stringIsProbablyANumber(effortString)) {
		// case 1: is a number: use it for all covariates
		for (size_t c = 0; c < numCovariatesEffort(); ++c) { _simulateOneCovariateEffort(c, effortString); }
	} else if (str::stringContains(effortString, ',')) {
		// cases 2 and 4: comma-separated string of numbers (case 2) or distributions (case 4) or both
		std::vector<std::string> vec = _splitVecWithDistributions(effortString);
		if (vec.size() != numCovariatesEffort()) {
			UERROR("The size of covariates provided with argument 'covariatesEffort' (", vec.size(),
			       ") does not match the number of covariates (", numCovariatesEffort(), ").");
		}
		for (size_t c = 0; c < vec.size(); ++c) { _simulateOneCovariateEffort(c, vec[c]); }
	} else {
		// case 3: use the same distribution for all covariates
		for (size_t c = 0; c < numCovariatesEffort(); ++c) { _simulateOneCovariateEffort(c, effortString); }
	}
}

void TMethods::_simulateOneCovariateDetection(size_t c, const std::string &Cov) {
	if (coretools::str::stringIsProbablyANumber(Cov)) {
		// case 1: is a number: use it for all covariates
		double value = coretools::str::fromString<double>(Cov);
		for (auto &_location : _locations) {
			for (auto &t : _location) { t.setCovariateDetection(c, value); }
		}
	} else {
		// case 2: parse distributions. "normal(0, 1), uniform(a,b)
		const auto [distr, params] = _getParamsAndDistr(Cov);
		if (distr == "normal") {
			_simulateCovariateFromDistribution<coretools::probdist::TNormalDistr, false>(c, params);
		} else if (distr == "uniform") {
			_simulateCovariateFromDistribution<coretools::probdist::TUniformDistr, false>(c, params);
		} else {
			UERROR("Unknown distribution '", distr,
			       "' (argument 'covariatesDetection'). Supported distributions are: normal and uniform.");
		}
	}
}

void TMethods::_simulateCovariatesDetection() {
	using namespace coretools;
	std::string effortString = instances::parameters().getParameterWithDefault("covariatesDetection", "normal(0, 1)");
	if (str::stringIsProbablyANumber(effortString)) {
		// case 1: is a number: use it for all covariates
		for (size_t c = 0; c < numCovariatesDetection(); ++c) { _simulateOneCovariateDetection(c, effortString); }
	} else if (str::stringContains(effortString, ',')) {
		// cases 2 and 4: comma-separated string of numbers (case 2) or distributions (case 4) or both
		std::vector<std::string> vec = _splitVecWithDistributions(effortString);
		if (vec.size() != numCovariatesDetection()) {
			UERROR("The size of covariates provided with argument 'covariatesDetection' (", vec.size(),
			       ") does not match the number of covariates (", numCovariatesDetection(), ").");
		}
		for (size_t c = 0; c < vec.size(); ++c) { _simulateOneCovariateDetection(c, vec[c]); }
	} else {
		// case 3: use the same distribution for all covariates
		for (size_t c = 0; c < numCovariatesDetection(); ++c) { _simulateOneCovariateDetection(c, effortString); }
	}
}

void TMethods::_simulateCovariates() {
	_simulateCovariatesEffort();
	_simulateCovariatesDetection();

	// standarize covariates
	standardizeCovariates();

	// now that we have simulated all covariates: re-calculate effort
	initializeEffort();
}

void TMethods::_setValuesForSingleCovariate() {
	if (!inferAlpha()) {
		for (size_t c = 0; c < _alpha->totalSize(); ++c) { _alpha->set(c, 1.0); }
		_alpha->setIsUpdated(false);
	}
	if (!inferBeta()) {
		_beta0->set(0.0); // set beta0 to zero
		_beta0->setIsUpdated(false);

		for (size_t c = 0; c < _beta->totalSize(); ++c) { _beta->set(c, 0.0); }
		_beta->setIsUpdated(false);
	}
}

void TMethods::simulatePoisson(size_t N0, const TModelBase &Model, const TUniqueContainer<std::string> &locations,
                               bool SimulateCovariates) {
	// simulate covariates, if necessary
	_alpha->normalizeToSum1AfterUpdate();
	_setValuesForSingleCovariate();

	if (SimulateCovariates) { _simulateCovariates(); }

	for (size_t j = 0; j < size(); j++) { // loop over all locations
		_locations[j].simulatePoisson(N0, Model);
	}
}

std::vector<double> TMethods::simulateNB(size_t N0, const TModelBase &Model,
                                         const TUniqueContainer<std::string> &locations, double A,
                                         bool SimulateCovariates) {
	// simulate covariates, if necessary
	_alpha->normalizeToSum1AfterUpdate();
	_setValuesForSingleCovariate();

	if (SimulateCovariates) { _simulateCovariates(); }

	std::vector<double> mus(_locations.size());
	for (size_t j = 0; j < size(); j++) { // loop over all locations
		mus[j] = (double)N0;
		_locations[j].simulateNB(Model, A, mus[j]);
	}
	return mus;
}

void TMethods::write(coretools::TOutputFile &File, const TUniqueContainer<std::string> &Locations,
                     const TUniqueContainer<size_t> &Timepoints) const {
	for (const auto &location : _locations) { location.write(File, Locations, Timepoints); }
}

void TMethods::estimateInitialBeta() {
	// TODO: How to initialize alpha and beta?
	// set initial value of all alpha's and beta's to one
	for (size_t c = 0; c < _alpha->totalSize(); ++c) { _alpha->setInitial(c, 1.0); }
	for (size_t c = 0; c < _beta->totalSize(); ++c) { _beta->setInitial(c, 1.0); }
	_beta0->setInitial(0.0);

	_alpha->normalizeToSum1AfterUpdate();

	// set effort based on initial values of beta
	initializeEffort();
	_setValuesForSingleCovariate();
}

void TMethods::_updateBeta(size_t c1, size_t c2, const TData &Data,
                           const std::vector<std::unique_ptr<TModelBase>> &Models) {
	// propose new value
	_beta->updatePairWiseSumOneVector(c1, c2);

	// update effort for all locations and timepoints of current method
	for (auto &location : _locations) {
		location.updateEffort(_beta->value(c1), _beta->oldValue(c1), _beta->value(c2), _beta->oldValue(c2), c1, c2);
	}
	// calculate Hastings ratio
	double logRatio = _beta->getLogPriorRatio(c1) + _beta->getLogPriorDensity(c2);
	for (size_t s = 0; s < Models.size(); ++s) { // loop over species
		logRatio += Models[s]->calculateLLRatio_perMethod(method_id(), Data);
	}
	// decide if accept or reject
	coretools::TRange range(c1, c2);
	if (_beta->acceptOrReject(logRatio, range, c1)) {
		// accepted: set likelihood in all epoch models
		for (size_t s = 0; s < Models.size(); ++s) { Models[s]->swapTryCur_perMethod(method_id(), Data); }
	} else {
		// rejected: re-set effort to old value before update
		for (auto &location : _locations) { location.resetEffort(); }
	}
}

void TMethods::updateBeta(const TData &Data, const std::vector<std::unique_ptr<TModelBase>> &Models) {
	if (!_beta->isUpdated()) { return; }

	stattools::TIndexSampler sampler(numCovariates()); // draw random indices from range [0, numCovariates]
	sampler.sampleIndices();

	for (size_t c = 0; c < sampler.length(); ++c) {
		// pick two indices i to update
		auto [c1, c2] = sampler.getIndexPair(c);
		_updateBeta(c1, c2, Data, Models);
	}
}

//-------------------------------------------
// TData
//-------------------------------------------

TData::TData() = default;

void TData::addMethod(const TMethods &Method) { _methods.emplace_back(Method); }

void TData::_fillNumMethLocs() {
	// check how many combinations of methods/locations we have - used to initialize vector of LL
	_numMethLoc.resize(_numSpecies, 0);
	for (size_t species = 0; species < _numSpecies; ++species) {
		for (const auto &m : _methods) {
			if (m.hasDataForSpeciesID(species)) { _numMethLoc[species] += m.size(); }
		}
	}
}

std::pair<size_t, size_t> TData::_getMethAndLocIndex(size_t ij_index, size_t speciesID) const {
	// returns the actual indices of specific method and location, not id in unique container.
	size_t counter = 0;
	for (size_t i = 0; i < _methods.size(); i++) {
		if (_methods[i].hasDataForSpeciesID(speciesID)) {
			if (counter + _methods[i].size() > ij_index) {
				size_t j = ij_index - counter;
				return {i, j};
			}
		}
		counter += _methods[i].size();
	}
	DEVERROR("Did not find method and location index for species - this should not happen.");
}

void TData::fillMethLocIndices(size_t NumSpecies) {
	_numSpecies = NumSpecies;

	_fillNumMethLocs();

	// fill linear to i_j
	_linear_to_i_j.resize(_numSpecies);
	for (size_t species = 0; species < _numSpecies; ++species) {
		_linear_to_i_j[species].resize(_numMethLoc[species]);
		for (size_t ij = 0; ij < _numMethLoc[species]; ++ij) {
			_linear_to_i_j[species][ij] = _getMethAndLocIndex(ij, species);
		}
	}

	// fill i_j to linear
	_i_j_to_linear.resize(_numSpecies);
	for (size_t species = 0; species < _numSpecies; ++species) {
		_i_j_to_linear[species].resize(_methods.size());           // size = numMethods
		for (size_t i = 0; i < _methods.size(); ++i) {
			_i_j_to_linear[species][i].resize(_methods[i].size()); // size = numLocations of this method
			for (size_t j = 0; j < _methods[i].size(); ++j) {
				auto it =
				    std::find(_linear_to_i_j[species].begin(), _linear_to_i_j[species].end(), std::make_pair(i, j));
				_i_j_to_linear[species][i][j] = std::distance(_linear_to_i_j[species].begin(), it);
			}
		}
	}
}

size_t TData::size() const { return _methods.size(); }

const TMethods &TData::operator[](const size_t index) const { return _methods[index]; }

TMethods &TData::operator[](const size_t index) { return _methods[index]; }

std::vector<TMethods>::const_iterator TData::cbegin() const { return _methods.cbegin(); }

std::vector<TMethods>::iterator TData::begin() { return _methods.begin(); }

std::vector<TMethods>::const_iterator TData::cend() const { return _methods.cend(); }

std::vector<TMethods>::iterator TData::end() { return _methods.end(); }
TMethods &TData::back() { return _methods.back(); }

size_t TData::numMethLoc(size_t SpeciesIDInUniqueContainer) const { return _numMethLoc[SpeciesIDInUniqueContainer]; }
size_t TData::i_j_to_linear(size_t SpeciesIDInUniqueContainer, size_t i, size_t j) const {
	return _i_j_to_linear[SpeciesIDInUniqueContainer][i][j];
}
const std::pair<size_t, size_t> &TData::linear_to_i_j(size_t SpeciesIDInUniqueContainer, size_t linear) const {
	return _linear_to_i_j[SpeciesIDInUniqueContainer][linear];
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