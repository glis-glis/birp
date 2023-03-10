//
// Created by liam on 10.03.23.
//

#ifndef BIRP_TDATA_H
#define BIRP_TDATA_H

#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "stattools/ParametersObservations/TParameterTyped.h"
#include <vector>
#include <string>


//-------------------------------------------
// TTimepoints
//-------------------------------------------

class TTimepoints {
private:

public:
    std::vector<double> _covariates;
    std::vector<size_t> _counts_per_species;
    TypeEffort effort   = 0;
    TypeEffort effort_old   = 0;
    size_t timepoint_id = 0;
    size_t _covariateName_id;

    TTimepoints() = default;

    // setters

    // getters
    size_t size() const;


};



//-------------------------------------------
// TLocations
//-------------------------------------------

class TLocations {
private:
    std::string _locationName;
    std::vector<TTimepoints> _timepoints;

public:
    size_t location_id = 0;
    size_t covariate_id = 0;
    TLocations() = default;

    // setters

    // getters
    size_t size() const;

    const TTimepoints &operator[](const size_t index) const { return _timepoints[index]; }

    std::vector<TTimepoints>::const_iterator cbegin() const;
    std::vector<TTimepoints>::iterator begin();
    std::vector<TTimepoints>::const_iterator cend() const;
    std::vector<TTimepoints>::iterator end();
};


//-------------------------------------------
// TMethods
//-------------------------------------------

class TMethods {
private:
    std::string _methodName;
    std::vector<std::string> _speciesNames;
    std::vector<double> _covariates;
    std::vector<TLocations> _locations;

public:
    explicit TMethods(const std::string &MethodsName);
    TMethods() = default;

    // setters

    // getters
    size_t size() const;

    const TLocations &operator[](const size_t index) const { return _locations[index]; }

    std::vector<TLocations>::const_iterator cbegin() const;
    std::vector<TLocations>::iterator begin();
    std::vector<TLocations>::const_iterator cend() const;
    std::vector<TLocations>::iterator end();
};

#endif //BIRP_TDATA_H
