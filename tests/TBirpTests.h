//
// Created by caduffm on 1/5/23.
//

#ifndef BIRP_TBIRPTESTS_H
#define BIRP_TBIRPTESTS_H

#include "TBirpCore.h"
#include "TData.h"

class TBirpCoreBridge : public TBirpCore {
	// Dummy class inheriting from TBirpCore: allows for testing protected functions
public:
	void _readFile(const std::string &Filename) { return TBirpCore::_readFile(Filename, true); }
	void _sortTimepoints() { return TBirpCore::_sortTimepoints(); }

	// getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string> &locations() const { return _locations; };
	const std::vector<TMethods> &methods() const { return _methods; };
	const TUniqueContainer<size_t> &timepoints() const { return _timepoints; };
    const TUniqueContainer<std::string> &speciesNames() const {return _speciesNames;}
};

#endif // BIRP_TBIRPTESTS_H
