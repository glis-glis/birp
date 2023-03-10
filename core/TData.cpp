//
// Created by liam on 10.03.23.
//

#include "TData.h"

//-------------------------------------------
// TTimepoints
//-------------------------------------------

size_t TTimepoints::timepoint_id(){return _timepoint_id;}

//-------------------------------------------
// TLocations
//-------------------------------------------

size_t TLocations::size() const{return _timepoints.size();}
size_t TLocations::location_id(){return _location_id;}

const &TLocations::operator[](const size_t index) const { return _timepoints[index]; }  // TODO ??

void TLocations::add(const TTimepoints &pt) {
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
std::vector<std::string> TMethods::speciesNames(){return _speciesNames;}

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

std::vector<TLocations>::const_iterator TMethods::cbegin() const { return _locations.cbegin(); }
std::vector<TLocations>::iterator TMethods::begin() { return _locations.begin(); }
std::vector<TLocations>::const_iterator TMethods::cend() const { return _locations.cend(); }
std::vector<TLocations>::iterator TMethods::end() { return _locations.end(); }