//
// Created by liam on 3/01/23.
//

#include "TSpecies.h"
#include "BirpTypes.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include "coretools/Storage/TStorage.h"
#include "coretools/Strings/splitters.h"
#include "coretools/algorithms.h"
#include <algorithm>
#include <vector>
#include "coretools/Main/TLog.h"

TDataPoint::TDataPoint(size_t Timpoint_id, TypeCounts Counts, TypeEffort Effort) {
	timepoint_id = Timpoint_id;
	counts       = Counts;
	effort       = Effort;
}

TDataPoint::TDataPoint(size_t Timpoint_id) {
	timepoint_id = Timpoint_id;
}

void TDataPoint::_simulatePoisson(double lambda){
	lambda *= effort;
	if (lambda == 0.0) {
		counts = 0;
	} else {
		counts = coretools::instances::randomGenerator().getPoissonRandom((double)lambda);
	}
}

void TDataPoint::simulate(size_t N0, double phi, double detectionProbs) {
	coretools::Positive<double> lambda = (double)N0 * phi * detectionProbs;
	_simulatePoisson(lambda);
}
void TDataPoint::simulate(double Effort, double phi) {
	effort = Effort;
	_simulatePoisson(phi);
}


bool TDataPoint::operator<(const TDataPoint &other) const { return timepoint_id < other.timepoint_id; }




TDataSeries::TDataSeries(size_t Location_id, size_t Method_id) {
	_location_id = Location_id;
	_method_id   = Method_id;
}

TDataSeries::TDataSeries(const TDataPoint &Datapoint, size_t Location_id, size_t Method_id) {
	_location_id = Location_id;
	_method_id   = Method_id;
	add(Datapoint);
}

TDataSeries::TDataSeries(size_t Location_id, size_t Method_id, size_t numTimepoints) {
	_location_id = Location_id;
	_method_id   = Method_id;
	_dataPoints.reserve(numTimepoints);
	for(size_t i = 0; i < numTimepoints; i++){
		_dataPoints.emplace_back(i);
	}
}

// set
void TDataSeries::add(const TDataPoint &Datapoint) {
	_dataPoints.push_back(Datapoint);
	_nu += Datapoint.counts;
}

// get
size_t TDataSeries::size() const { return _dataPoints.size(); }


std::vector<TDataPoint>::const_iterator TDataSeries::cbegin() const { return _dataPoints.cbegin(); }
std::vector<TDataPoint>::iterator TDataSeries::begin() { return _dataPoints.begin(); }
std::vector<TDataPoint>::const_iterator TDataSeries::cend() const { return _dataPoints.cend(); }
std::vector<TDataPoint>::iterator TDataSeries::end() { return _dataPoints.end(); }

// getters
size_t TDataSeries::location_id() const { return _location_id; }
size_t TDataSeries::method_id() const { return _method_id; }
TypeCounts TDataSeries::nu() const { return _nu; }

void TDataSeries::sorttimes(const std::vector<size_t> &index) {
	std::vector<size_t> positionOfTimepointIds;
	for (size_t i = 0; i < _dataPoints.size(); i++) {
		auto it    = std::find(index.begin(), index.end(), _dataPoints[i].timepoint_id);
		size_t pos = it - index.begin();
		positionOfTimepointIds.push_back(pos);
	}
	// sort the newly found positions of timepoints using ranksort()
	std::vector<size_t> sortingIndex = coretools::rankSort(positionOfTimepointIds);
	// sort _dataPoints
	_dataPoints                      = coretools::sortContainerByRank(_dataPoints, sortingIndex);
	// now sort positionOfTimepointIds
	positionOfTimepointIds           = coretools::sortContainerByRank(positionOfTimepointIds, sortingIndex);
	// now set timepoint_id
	for (size_t i = 0; i < _dataPoints.size(); i++) { _dataPoints[i].timepoint_id = positionOfTimepointIds[i]; }
}

void TDataSeries::checkData(const std::string &name, const TUniqueContainer<std::string> &locations,
                            const TUniqueContainer<std::string> &methods) {
	size_t nonZeroEffort = 0;
	for (size_t i = 0; i < _dataPoints.size(); i++) {
		if (_dataPoints[i].effort > 0) { nonZeroEffort += 1; }
	}
	if (nonZeroEffort < 2) {
		coretools::instances::logfile().warning(
		    "The location ", locations[_location_id], " with method ", methods[_method_id], " for species ", name,
		    " contains less than two timepoints with non-zero effort, therefore it won't be analysed.");
	}
}

void TDataSeries::simulate(size_t N0, const std::vector<double>& phi, double detectionProbs) {
	for (size_t t = 0; t < _dataPoints.size(); t++) { _dataPoints[t].simulate(N0, phi[t], detectionProbs); }
}
void TDataSeries::simulate(double effort, const std::vector<double>& phi){
	for (size_t t = 0; t < _dataPoints.size(); t++) { _dataPoints[t].simulate(effort, phi[t]); }
}



TSpecies::TSpecies(const std::string &Name) { _name = Name; }
TSpecies::TSpecies(const std::string &Name, const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex) {
	_name = Name;
	add(pt, LocationIndex, MethodIndex);
}
TSpecies::TSpecies(const std::string &Name, size_t numLocation, size_t numMethods, size_t numTimepoints){
	_name = Name;
	_data.reserve(numLocation * numMethods);
	for(size_t i = 0; i < numLocation; i++){
		for(size_t j = 0; j < numMethods; j++){
			_data.emplace_back(i,j, numTimepoints);
		}
	}
}

std::string TSpecies::name() const { return _name; }

void TSpecies::add(const TDataPoint &pt, size_t LocationIndex, size_t MethodIndex) {
	bool dataPointExists = false;
	for (size_t d = 0; d < _data.size(); d++) {
		if (_data[d].location_id() == LocationIndex && _data[d].method_id() == MethodIndex) {
			dataPointExists = true;
			_data[d].add(pt);
			break;
		}
	}
	if (!dataPointExists) { _data.emplace_back(pt, LocationIndex, MethodIndex); }
}

void TSpecies::sorttimes(const std::vector<size_t> &index) {
	for (size_t ds = 0; ds < _data.size(); ds++) { _data[ds].sorttimes(index); }
}

size_t TSpecies::size() const { return _data.size(); }

void TSpecies::checkData(const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods) {
	for (size_t i = 0; i < _data.size(); i++) { _data[i].checkData(_name, locations, methods); }
}

std::vector<TDataSeries>::const_iterator TSpecies::cbegin() const { return _data.cbegin(); }
std::vector<TDataSeries>::iterator TSpecies::begin() { return _data.begin(); }
std::vector<TDataSeries>::const_iterator TSpecies::cend() const { return _data.cend(); }
std::vector<TDataSeries>::iterator TSpecies::end() { return _data.end(); }


void TSpecies::simulate(size_t N0, const std::vector<double> &phi,
                        const coretools::TMultiDimensionalStorage<double, 2> &detectionProbs,
                        const TUniqueContainer<std::string> &locations, const TUniqueContainer<std::string> &methods) {
	for (size_t i = 0; i < size(); i++) {
		std::string locName  = locations[_data[i].location_id()];
		std::string methName = methods[_data[i].method_id()];
		size_t locIndex      = detectionProbs.getDimensionName(0)->getIndex(locName);
		size_t methIndex     = detectionProbs.getDimensionName(1)->getIndex(methName);
		_data[i].simulate(N0, phi, detectionProbs(locIndex, methIndex));
	}
}

void TSpecies::simulate(const std::vector<double>& phi) {
	std::string effortString = coretools::instances::parameters().getParameterWithDefault("effort", "Poisson(10)");
	if (coretools::str::stringIsProbablyANumber(effortString)) {
		// case 1: is a number: use it for all efforts
		double effort = coretools::str::fromString<double>(effortString);
		for (size_t i = 0; i < size(); i++) {
			_data[i].simulate(effort, phi);
		}
	} else {
		// case 2: parse distributions. "Poisson(1), uniform(a,b), zero-poisson(p, lambda) -> p: p(not 0),
		coretools::str::TSplitter sp(effortString, '(');
		std::string_view distr = sp.front();          // takes first element
		sp.popFront();                                // takes 'Poisson' away
		std::string params = (std::string)sp.front(); // now parameters are in front
		coretools::str::extractBefore(params, ')');
		if (distr == "Poisson") {
			coretools::probdist::TPoissonDistr poisDistr(params);
			for (size_t i = 0; i < size(); i++) {
				_data[i].simulate(poisDistr, phi);
			}
		} else if (distr == "uniform") {
			coretools::probdist::TUniformDistr uniformDist(params);
			for (size_t i = 0; i < size(); i++) {
				_data[i].simulate(uniformDist, phi);
			}
		} else if (distr == "zeroPoisson") {
			TZeroPoisDistr zeroPoisDist(params);
			for (size_t i = 0; i < size(); i++) {
				_data[i].simulate(zeroPoisDist, phi);
			}
		}
	}
}





//------------------------------------------------
// zero-Poisson distribution
//------------------------------------------------

TZeroPoisDistr::TZeroPoisDistr(coretools::Probability p, coretools::StrictlyPositive<double> lambda)
    : _p(p), _lambda(lambda) {
}

void TZeroPoisDistr::set(coretools::Probability p, coretools::StrictlyPositive<double> Lambda) {
	_lambda = Lambda;
	_p = p;
}

void TZeroPoisDistr::set(std::string_view parameterString) {
	coretools::str::convertString(parameterString, std::string{"Use "}.append(name) + "(p, lambda) with 0 <= p <= 1, lambda>0.", _p, _lambda);
}

static size_t sample(coretools::Probability p, coretools::StrictlyPositive<double> lambda) {
	double tmp = coretools::instances::randomGenerator().getRand(0, 1);
	if (tmp > p) { // p is probability of taking Poisson
		return 0;
	} else {
		return coretools::instances::randomGenerator().getPoissonRandom(lambda);
	}
}

size_t TZeroPoisDistr::sample() const noexcept { return sample(_p, _lambda); }
