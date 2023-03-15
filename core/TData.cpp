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

TLocations::TLocations(const TTimepoints &pt, size_t location_id){
    _location_id = location_id;
    add(pt);
}

const TTimepoints &TLocations::operator[](const size_t index) const { return _timepoints[index]; }

void TLocations::add(TTimepoints pt) {
    for (auto &t : _timepoints) {
        // TODO if timepoint already exists; throw error. same method, loc and time should not occur
        if (t.timepoint_id() == pt.timepoint_id()) {
            UERROR("Make sure your ipnut file does not contain multiple rows with identical method, location and timepoints");
        }
    }
    // when we are here: timpoint is uniqe -> push it back
    _timepoints.push_back(pt);
    // now add nu
    //TODO make loop over all species
    for(auto c:pt.counts_per_species()){
        _nu[c] += (uint32_t)pt.counts_per_species()[c];
    }

}

std::vector<TTimepoints>::const_iterator TLocations::cbegin() const { return _timepoints.cbegin(); }
std::vector<TTimepoints>::iterator TLocations::begin() { return _timepoints.begin(); }
std::vector<TTimepoints>::const_iterator TLocations::cend() const { return _timepoints.cend(); }
std::vector<TTimepoints>::iterator TLocations::end() { return _timepoints.end(); }

//-------------------------------------------
// TMethods
//-------------------------------------------

TMethods::TMethods(const std::string_view &methodName) { _methodName = methodName; }


size_t TMethods::size() const{ return _locations.size();}


void TMethods::add(const TTimepoints &pt, size_t locationIndex) {
    bool dataPointExists = false;
    for (auto &loc : _locations) {
        if (loc.location_id() == locationIndex) {
            dataPointExists = true;
            loc.add(pt);
            break;
        }
    }
    if (!dataPointExists) { _locations.emplace_back(pt, locationIndex); }
}

std::vector<TLocations>::const_iterator TMethods::cbegin() const { return _locations.cbegin(); }
std::vector<TLocations>::iterator TMethods::begin() { return _locations.begin(); }
std::vector<TLocations>::const_iterator TMethods::cend() const { return _locations.cend(); }
std::vector<TLocations>::iterator TMethods::end() { return _locations.end(); }