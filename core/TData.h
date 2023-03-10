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
	// TODO: mach doch die aui private und schrib getter/setter
    std::vector<double> _covariates;
    std::vector<size_t> _counts_per_species;
    TypeEffort effort   = 0;
    TypeEffort effort_old   = 0;
    size_t timepoint_id = 0;
    size_t _covariateName_id; // TODO: ds wär ä vektor: eine pro Element i de covariates

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
    std::string _locationName; // TODO: entweder name oder id, aber beides macht nid Sinn - mach doch eifach t'id
    std::vector<TTimepoints> _timepoints;

public:
    size_t location_id = 0;
    size_t covariate_id = 0; // TODO: dä bruchts hie nid
    TLocations() = default;

    // setters

    // getters
    size_t size() const;

    const TTimepoints &operator[](const size_t index) const { return _timepoints[index]; } // TODO: schieb aui implementatione iz cpp file, hie nume definitione

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
    std::vector<double> _covariates; // TODO: ghöre nid hie härä, mit "covariate coefficients" hani der parameter (t'beta) gmeint
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
