//
// Created by caduffm on 7/20/22.
//

#ifndef TBIRPPRIOR_H
#define TBIRPPRIOR_H

#include "BirpTypes.h"
#include "TMatrix.h"
#include "TObservationTyped.h"
#include "TParameterTyped.h"
#include "TSpecies.h"

template<typename T>
class TUniqueContainer{
private:
	std::vector<T> _names;
public:
	bool exists(const T &Value) const {
		return std::find(_names.begin(), _names.end(),Value)!=_names.end();
	}
	void add(const T &Value){
		if (!exists(Value)) {
			_names.push_back(Value);
		}
	}
	size_t getIndex(const T &Value) const {
		auto it = std::find(_names.begin(), _names.end(), Value);
		if (it != _names.end()) {
			return  it - _names.begin();
		}
		DEVERROR("Value ",Value, " not found in location-vector. This should not happen.");
	}

	size_t size() const {
		return _names.size();
	}

	const T &operator[](size_t Index) const { return _names[Index]; }
	typename std::vector<T>::const_iterator cbegin() const { return _names.cbegin(); }
	typename std::vector<T>::const_iterator cend() const { return _names.cend(); }

};


class TBirpPrior : public stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2> {
private:
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::Storage;
	using typename stattools::prior::TBaseLikelihoodPrior<stattools::TObservationBase, TypeCounts, 2>::UpdatedStorage;

	// Birp parameters and co-variates
	stattools::TParameterTyped<TypeGamma, 1> *_gamma     = nullptr;
	stattools::TObservationTyped<TypeEffort, 2> *_effort = nullptr;

	// main dimensions
	size_t _numMethods    = 0;
	size_t _numTimePoints = 0;

	// epochs, rho
	size_t _numEpochs = 0;
	coretools::TMatrix<size_t> _rho;
	std::vector<TypeCounts> _nu; // sums of counts at one location and method (ij) across all time points (k)

	// names
	TUniqueContainer<std::string> _locations;
	TUniqueContainer<std::string> _methods;
	TUniqueContainer<size_t> _timepoints;

	// temporary variables, for speed
	double _oldLL           = 0.0;
	double _oldLogJeffPrior = 0.0;

	// other stuff
	std::string _filename;

protected:
	// simulate values
	void _simulateUnderPrior(Storage *Data) override;

	// likelihood
	double _calculateLL(const coretools::TMatrix<double> & p) const;
	double _calculateLL() const;

	// Jeffrey's prior
	double _calculateLogJeffreysPrior(const coretools::TMatrix<double> & p) const;

	// MCMC updates
	void _evaluateUpdateGamma(size_t m);

public:
	TBirpPrior(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TObservationTyped<TypeEffort, 2> *Effort,
	           const std::string &Filename);
	~TBirpPrior() override = default;

	[[nodiscard]] std::string name() const override;

	// initialization
	void initializeInferred() override;

	// estimate initial values
	void estimateInitialPriorParameters() override;

	// full log densities
	double getSumLogPriorDensity(const Storage &) const override;

	// update all parameters
	void updateParams() override;

	// getters
	double rho(size_t Time, size_t Epoch) const;
	coretools::TMatrix<size_t> rhoMat() const;
};





/*// not needed anymore
class TLocations {
private:
	std::vector<std::string> _names;
	size_t _size;
public:
};

class TTimePoints {
private:
	std::vector<size_t> _timepoints;
	size_t _size;
public:
};

class TMethods {
private:
	std::vector<std::string> _names;
	size_t _size;
public:
};*/


#endif // TBIRPPRIOR_H
