//
// Created by madleina on 03.03.23.
//

#include "TCovariates.h"

//-------------------------------------------
// TCovariate_icjk
//-------------------------------------------

TCovariate_icjk::TCovariate_icjk(size_t Timpoint_id, TypeCovariate Covariate) {
	timepoint_id = Timpoint_id;
	covariate    = Covariate;
}

//-------------------------------------------
// TCovariate_icj
//-------------------------------------------

TCovariate_icj::TCovariate_icj(const TCovariate_icjk &Covariate_icjk, size_t Location_id) {
	_location_id = Location_id;
	addTimePoint(Covariate_icjk);
}

void TCovariate_icj::addTimePoint(const TCovariate_icjk &Covariate_icjk) { _covariates_icjk.push_back(Covariate_icjk); }

size_t TCovariate_icj::size() const { return _covariates_icjk.size(); }

std::vector<TCovariate_icjk>::const_iterator TCovariate_icj::cbegin() const { return _covariates_icjk.cbegin(); }

std::vector<TCovariate_icjk>::iterator TCovariate_icj::begin() { return _covariates_icjk.begin(); }

std::vector<TCovariate_icjk>::const_iterator TCovariate_icj::cend() const { return _covariates_icjk.cend(); }

std::vector<TCovariate_icjk>::iterator TCovariate_icj::end() { return _covariates_icjk.end(); }

size_t TCovariate_icj::location_id() const { return _location_id; }

void TCovariate_icj::checkCovExists(size_t timepoint_id) {
	bool timepointExists = false;
	for (auto c : _covariates_icjk) {
		if (timepoint_id == c.timepoint_id) {
			timepointExists = true;
			break;
		}
	}
	if (!timepointExists) {
		UERROR("Not all timepoints used in counts file are also present in effort (covariate) file.");
	}
}

//-------------------------------------------
// TCovariate_ic
//-------------------------------------------

TCovariate_ic::TCovariate_ic(const TCovariate_icjk &Covariate_icjk, size_t Covariate_id, size_t Location_id) {
	_covariate_id = Covariate_id;
	addTimePoint(Covariate_icjk, Location_id);
}

void TCovariate_ic::addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Location_id) {
	// check if location already exists: if yes, add timepoint, else create new location and add timepoint
	bool locationExists = false;
	for (auto &d : _covariates_icj) {
		if (d.location_id() == Location_id) {
			locationExists = true;
			d.addTimePoint(Covariate_icjk);
			break;
		}
	}
	if (!locationExists) { _covariates_icj.emplace_back(Covariate_icjk, Location_id); }
}

size_t TCovariate_ic::size() const { return _covariates_icj.size(); }

std::vector<TCovariate_icj>::const_iterator TCovariate_ic::cbegin() const { return _covariates_icj.cbegin(); }

std::vector<TCovariate_icj>::iterator TCovariate_ic::begin() { return _covariates_icj.begin(); }

std::vector<TCovariate_icj>::const_iterator TCovariate_ic::cend() const { return _covariates_icj.cend(); }

std::vector<TCovariate_icj>::iterator TCovariate_ic::end() { return _covariates_icj.end(); }

size_t TCovariate_ic::covariate_id() const { return _covariate_id; }

void TCovariate_ic::checkCovExists(size_t location_id, size_t timepoint_id) {
	bool locationExists = false;
	for (auto c : _covariates_icj) {
		if (location_id == c.location_id()) {
			locationExists = true;
			c.checkCovExists(timepoint_id);
			break;
		}
	}
	if (!locationExists) {
		UERROR("Not all locations used in counts file are also present in effort (covariate) file.");
	}
}

//-------------------------------------------
// TCovariate_i
//-------------------------------------------

TCovariate_i::TCovariate_i(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id,
                           size_t Location_id) {
	_method_id = Method_id;
	addTimePoint(Covariate_icjk, Covariate_id, Location_id);
}

void TCovariate_i::addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Covariate_id, size_t Location_id) {
	// check if location already exists: if yes, add timepoint, else create new location and add timepoint
	bool covariateExists = false;
	for (auto &d : _covariates_ic) {
		if (d.covariate_id() == Covariate_id) {
			covariateExists = true;
			d.addTimePoint(Covariate_icjk, Location_id);
			break;
		}
	}
	if (!covariateExists) { _covariates_ic.emplace_back(Covariate_icjk, Covariate_id, Location_id); }
}

void TCovariate_i::checkCovExists(size_t location_id, size_t timepoint_id) {
	for (auto c : _covariates_ic) { c.checkCovExists(location_id, timepoint_id); }
}

size_t TCovariate_i::size() const { return _covariates_ic.size(); }

bool TCovariate_i::covariateIsEffort() const {
	// for one method: if a single covariate -> this is simply the effort
	// else: learn regression coefficients of covariates to get effort
	return _covariates_ic.size() == 1;
}

std::vector<TCovariate_ic>::const_iterator TCovariate_i::cbegin() const { return _covariates_ic.cbegin(); }

std::vector<TCovariate_ic>::iterator TCovariate_i::begin() { return _covariates_ic.begin(); }

std::vector<TCovariate_ic>::const_iterator TCovariate_i::cend() const { return _covariates_ic.cend(); }

std::vector<TCovariate_ic>::iterator TCovariate_i::end() { return _covariates_ic.end(); }

size_t TCovariate_i::method_id() const { return _method_id; }

//-------------------------------------------
// TCovariate
//-------------------------------------------

TCovariate::TCovariate(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id,
                       size_t Location_id) {
	addTimePoint(Covariate_icjk, Method_id, Covariate_id, Location_id);
}

void TCovariate::addTimePoint(const TCovariate_icjk &Covariate_icjk, size_t Method_id, size_t Covariate_id,
                              size_t Location_id) {
	// check if location already exists: if yes, add timepoint, else create new location and add timepoint
	bool methodExists = false;
	for (auto &d : _covariates_i) {
		if (d.method_id() == Method_id) {
			methodExists = true;
			d.addTimePoint(Covariate_icjk, Covariate_id, Location_id);
			break;
		}
	}
	if (!methodExists) { _covariates_i.emplace_back(Covariate_icjk, Method_id, Covariate_id, Location_id); }
}

void TCovariate::IDexists(size_t method_id, size_t location_id, size_t timepoint_id) {
	bool methodExists = false;
	for (auto &c : _covariates_i) {
		if (c.method_id() == method_id) {
			methodExists = true;
			c.checkCovExists(location_id, timepoint_id);
			break;
		}
	}
	if (!methodExists) { UERROR("Not all methods used in counts file are also present in effort (covariate) file."); }
}

size_t TCovariate::size() const { return _covariates_i.size(); }

std::vector<TCovariate_i>::const_iterator TCovariate::cbegin() const { return _covariates_i.cbegin(); }

std::vector<TCovariate_i>::iterator TCovariate::begin() { return _covariates_i.begin(); }

std::vector<TCovariate_i>::const_iterator TCovariate::cend() const { return _covariates_i.cend(); }

std::vector<TCovariate_i>::iterator TCovariate::end() { return _covariates_i.end(); }