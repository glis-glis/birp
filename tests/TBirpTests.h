//
// Created by caduffm on 1/5/23.
//

#ifndef BIRP_TBIRPTESTS_H
#define BIRP_TBIRPTESTS_H

#include "TBirpCore.h"

class TBirpCoreBridge : public TBirpCore {
	// Dummy class inheriting from TBirpCore: allows for testing protected functions
public:
	void _readFile(const std::string &Filename) { return TBirpCore::_readFile(Filename, true); }
	void _sortTimepoints() { return TBirpCore::_sortTimepoints(); }

	// getters for locations, methods, timepoints, species for tests
	const TUniqueContainer<std::string> &locations() const { return _locations; };
	const TUniqueContainer<std::string> &methods() const { return _methods; };
	const TUniqueContainer<size_t> &timepoints() const { return _timepoints; };
	const std::vector<TSpecies> &species() const { return _species; };
};

#endif // BIRP_TBIRPTESTS_H
