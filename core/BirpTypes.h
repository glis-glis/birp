//
// Created by caduffm on 7/20/22.
//

#ifndef BIRPTYPES_H
#define BIRPTYPES_H

#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include <cstdint>

//-------------------------------
// typedefs
//-------------------------------

// types of parameters/observations
typedef coretools::UnsignedInt<uint16_t> TypeCounts;
typedef coretools::UnsignedInt<uint32_t> TypeNu;
typedef coretools::UnsignedInt<uint16_t> TypeIndex;
typedef coretools::WeakType<double> TypeGamma;
typedef coretools::Positive<double> TypeEffort;
typedef coretools::UnsignedInt<uint16_t> TypeTime;
typedef coretools::StrictlyPositiveWithMax<double> TypeMu;
typedef coretools::StrictlyPositive<double> TypeSigmaMu;
typedef coretools::MinMaxVariable<double> TypeLogA;

#endif // BIRPTYPES_H
