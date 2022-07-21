//
// Created by caduffm on 7/20/22.
//

#ifndef BIRPTYPES_H
#define BIRPTYPES_H

#include "commonWeakTypes.h"
#include "probability.h"
#include <cstdint>

//-------------------------------
// typedefs
//-------------------------------

// types of parameters/observations
typedef coretools::UnsignedInt<uint16_t> TypeCounts;
typedef coretools::WeakType<double> TypeGamma;
typedef coretools::Positive<double> TypeEffort;

#endif // BIRPTYPES_H
