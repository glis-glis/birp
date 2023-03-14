//
// Created by liam on 10.03.23.
//

#include "TData.h"

//-------------------------------------------
// TTimepoints
//-------------------------------------------

TTimepoints::TTimepoints(std::vector<TypeCounts> countVec, std::vector<TypeCovariate> covVec){
    _covariates = covVec;
    _counts_per_species = countVec;
}

TTimepoints::TTimepoints(size_t timepoint_id){
    _timepoint_id = timepoint_id;
}

std::vector<TypeCovariate> TTimepoints::covariates(){return _covariates;}
std::vector<size_t> TTimepoints::covariate_id(){return _covariate_id;}
std::vector<TypeCounts> TTimepoints::counts_per_species(){return _counts_per_species;}
std::vector<size_t> TTimepoints::species_id(){return _species_id;}
TypeEffort TTimepoints::effort(){return _effort;}
TypeEffort TTimepoints::effort_old(){return _effort_old;}
size_t TTimepoints::timepoint_id(){return _timepoint_id;}


//-------------------------------------------
// TLocations
//-------------------------------------------

size_t TLocations::size() const{return _timepoints.size();}
size_t TLocations::location_id(){return _location_id;}

const &TLocations::operator[](const size_t index) const { return _timepoints[index]; }  // TODO ??

void TLocations::add(const TTimepoints &pt, size_t timepoint_id) {
    bool dataPointExists = false;
    for (auto &t : _timepoints) {
        if (t.timepoint_id() == timepoint_id) {
            dataPointExists = true;
            t.push_back(pt);
            _nu += (uint32_t)t.counts_per_species()[0]; // TODO which species to add to?
            break;
        }
    }
    if (!dataPointExists) { _timepoints.emplace_back(pt, timepoint_id); }  // TODO how to add constructor


    _timepoints.push_back(pt);
    _nu += (uint32_t)pt.counts; //TODO how to sum up per species?
}

std::vector<TTimepoints>::const_iterator TLocations::cbegin() const { return _timepoints.cbegin(); }
std::vector<TTimepoints>::iterator TLocations::begin() { return _timepoints.begin(); }
std::vector<TTimepoints>::const_iterator TLocations::cend() const { return _timepoints.cend(); }
std::vector<TTimepoints>::iterator TLocations::end() { return _timepoints.end(); }

//-------------------------------------------
// TMethods
//-------------------------------------------

TMethods::TMethods(const std::string &MethodsName) { _methodName = MethodsName; }


size_t TMethods::size() const{ return _locations.size();}


void TMethods::add(const TTimepoints &pt, size_t locationIndex, size_t timepointIndex) {  //TODO location and timepoints are not indices?
    bool dataPointExists = false;
    for (auto &loc : _locations) {
        if (loc.location_id() == locationIndex) {
            dataPointExists = true;
            loc.add(pt, timepointIndex);
            break;
        }
    }
    if (!dataPointExists) { _locations.emplace_back(pt, location_id); } //TODO  write constructor
}

std::vector<TLocations>::const_iterator TMethods::cbegin() const { return _locations.cbegin(); }
std::vector<TLocations>::iterator TMethods::begin() { return _locations.begin(); }
std::vector<TLocations>::const_iterator TMethods::cend() const { return _locations.cend(); }
std::vector<TLocations>::iterator TMethods::end() { return _locations.end(); }