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
    std::vector<TypeCovariate> _covariates;
    std::vector<size_t> _covariate_id;
    std::vector<TypeCounts> _counts_per_species;
    std::vector<size_t> _species_id;
    TypeEffort _effort   = 0;
    TypeEffort _effort_old   = 0;
    size_t _timepoint_id = 0;



public:
    TTimepoints() = default;
    TTimepoints(std::vector<TypeCounts> countVec, std::vector<TypeCovariate> covVec);
    TTimepoints(size_t timepoint_id);

    // setters

    // getters
    std::vector<TypeCovariate> covariates();
    std::vector<size_t> covariate_id();
    std::vector<TypeCounts> counts_per_species();
    std::vector<size_t> species_id();
    TypeEffort effort();
    TypeEffort effort_old();
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
    void add(const TTimepoints &pt, size_t timepoint_id);

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
    size_t _method_id;
    std::vector<double> _beta;
    std::vector<TLocations> _locations;

public:
    explicit TMethods(const std::string &MethodsName);
    TMethods() = default;

    // setters

    // getters
    size_t size() const;

    void add(const TTimepoints &pt, size_t location_id, size_t timepoint_id);

    const TLocations &operator[](const size_t index) const { return _locations[index]; }

    std::vector<TLocations>::const_iterator cbegin() const;
    std::vector<TLocations>::iterator begin();
    std::vector<TLocations>::const_iterator cend() const;
    std::vector<TLocations>::iterator end();
};

#endif //BIRP_TDATA_H
