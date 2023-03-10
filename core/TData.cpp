//
// Created by liam on 10.03.23.
//

#include "TData.h"

//-------------------------------------------
// TTimepoints
//-------------------------------------------



//-------------------------------------------
// TLocations
//-------------------------------------------

std::vector<TTimepoints>::const_iterator TLocations::cbegin() const { return _timepoints.cbegin(); }
std::vector<TTimepoints>::iterator TLocations::begin() { return _timepoints.begin(); }
std::vector<TTimepoints>::const_iterator TLocations::cend() const { return _timepoints.cend(); }
std::vector<TTimepoints>::iterator TLocations::end() { return _timepoints.end(); }

//-------------------------------------------
// TMethods
//-------------------------------------------

TMethods::TMethods(const std::string &MethodsName) { _methodName = MethodsName; }

size_t TMethods::size() const{ return _locations.size();}

std::vector<TLocations>::const_iterator TMethods::cbegin() const { return _locations.cbegin(); }
std::vector<TLocations>::iterator TMethods::begin() { return _locations.begin(); }
std::vector<TLocations>::const_iterator TMethods::cend() const { return _locations.cend(); }
std::vector<TLocations>::iterator TMethods::end() { return _locations.end(); }