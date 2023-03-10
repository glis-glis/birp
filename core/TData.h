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
    std::vector<double> _covariates;
    std::vector<size_t> _counts_per_species;
    TypeEffort _effort   = 0;
    TypeEffort _effort_old   = 0;
    size_t _timepoint_id = 0;
    std::vector<size_t> _covariateName_id; // TODO: ds wär ä vektor: eine pro Element i de covariates. so?

public:
    TTimepoints() = default;

    // setters

    // getters
    size_t timepoint_id();

};



//-------------------------------------------
// TLocations
//-------------------------------------------

class TLocations {
private:
    size_t _location_id = 0;
    std::vector<TTimepoints> _timepoints;
    std::vector<TypeNu> _nu;

public:
    TLocations() = default;

    // setters
    void add(const TTimepoints &pt);

    // getters
    size_t size() const;
    size_t location_id();

    const TLocations &operator[](const size_t index) const;


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
    std::vector<double> _beta;
    std::vector<TLocations> _locations;
    // TODO species_id?

public:
    explicit TMethods(const std::string &MethodsName);
    TMethods() = default;

    // setters

    // getters
    size_t size() const;
    std::vector<std::string> speciesNames();

    void add(const TTimepoints &pt, size_t location_id);

    const TLocations &operator[](const size_t index) const { return _locations[index]; }

    std::vector<TLocations>::const_iterator cbegin() const;
    std::vector<TLocations>::iterator begin();
    std::vector<TLocations>::const_iterator cend() const;
    std::vector<TLocations>::iterator end();
};

#endif //BIRP_TDATA_H
