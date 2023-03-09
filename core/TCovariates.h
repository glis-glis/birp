//
// Created by madleina on 03.03.23.
//

#ifndef BIRP_TCOVARIATES_H
#define BIRP_TCOVARIATES_H

#include "BirpTypes.h"
#include "TUniqueContainer.h"
#include "stattools/ParametersObservations/TParameterTyped.h"

//-------------------------------------------
// TCovariate_icjk
//-------------------------------------------

class TCovariate_icjk {
	// this class stores the value of method i and covariate c at location j at timepoint k
public:
	size_t timepoint_id     = 0;
	TypeCovariate covariate = 0.0;

	TCovariate_icjk(size_t Timpoint_id, TypeCovariate Covariate);
	TCovariate_icjk() = default;
};

//-------------------------------------------
// TCovariate_icj
//-------------------------------------------

class TCovariate_icj {
	// this class stores a vector over all timepoints for method i and covariate c at location j
private:
	std::vector<TCovariate_icjk> _covariates_icjk;
	size_t _location_id = 0;

public:
	explicit TCovariate_icj(const TCovariate_icjk &Covariate_icjk, size_t Location_id);

	// setters
	void addTimePoint(const TCovariate_icjk &Covariate_icjk);

	void checkCovExists(size_t timepoint);

	// getters
	size_t size() const;
	size_t location_id() const;

	const TCovariate_icjk &operator[](const size_t index) const { return _covariates_icjk[index]; }

	std::vector<TCovariate_icjk>::const_iterator cbegin() const;
	std::vector<TCovariate_icjk>::iterator begin();
	std::vector<TCovariate_icjk>::const_iterator cend() const;
	std::vector<TCovariate_icjk>::iterator end();
};

//-------------------------------------------
// TCovariate_ic
//-------------------------------------------

class TCovariate_ic {
	// this class stores a vector over all locations for method i and covariate c
private:
	std::vector<TCovariate_icj> _covariates_icj;
	size_t _covariate_id = 0;

public:
	explicit TCovariate_ic(const TCovariate_icjk &Covariate_icjk, size_t Covariate_id, size_t Location_id);

	// setters
	void addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Location_id);

	void checkCovExists(size_t location_id, size_t timepoint_id);

	// getters
	size_t size() const;
	size_t covariate_id() const;

	const TCovariate_icj &operator[](const size_t index) const { return _covariates_icj[index]; }

	std::vector<TCovariate_icj>::const_iterator cbegin() const;
	std::vector<TCovariate_icj>::iterator begin();
	std::vector<TCovariate_icj>::const_iterator cend() const;
	std::vector<TCovariate_icj>::iterator end();
};

//-------------------------------------------
// TCovariate_i
//-------------------------------------------

class TCovariate_i {
	// this class stores a vector over all covariates for method i
private:
	std::vector<TCovariate_ic> _covariates_ic;
	size_t _method_id       = 0;

public:
	explicit TCovariate_i(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id,
	                      size_t Location_id);

	// setters
	void addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Covariate_id, size_t Location_id);

	// getters
	size_t size() const;
	size_t method_id() const;
	bool covariateIsEffort() const;
	void checkCovExists(size_t location_id, size_t timepoint_id);

	const TCovariate_ic &operator[](const size_t index) const { return _covariates_ic[index]; }

	std::vector<TCovariate_ic>::const_iterator cbegin() const;
	std::vector<TCovariate_ic>::iterator begin();
	std::vector<TCovariate_ic>::const_iterator cend() const;
	std::vector<TCovariate_ic>::iterator end();
};

//-------------------------------------------
// TCovariate
//-------------------------------------------

class TCovariate {
	// this class stores a vector over all methods
private:
	std::vector<TCovariate_i> _covariates_i;

public:
	explicit TCovariate(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id,
	                    size_t Location_id);
	TCovariate() = default;

	// setters
	void addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id, size_t Location_id);

	// checkers
	void IDexists(size_t method_id, size_t location_id, size_t timepoint_id);

	// getters
	size_t size() const;

	const TCovariate_i &operator[](const size_t index) const { return _covariates_i[index]; }

	std::vector<TCovariate_i>::const_iterator cbegin() const;
	std::vector<TCovariate_i>::iterator begin();
	std::vector<TCovariate_i>::const_iterator cend() const;
	std::vector<TCovariate_i>::iterator end();
};

#endif // BIRP_TCOVARIATES_H
