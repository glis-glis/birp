//
// Created by caduffm on 7/20/22.
//

#ifndef BIRPTYPES_H
#define BIRPTYPES_H

#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include "coretools/Types/weakTypesWithLogExp.h"
#include <cstdint>

//-------------------------------
// typedefs
//-------------------------------

// types of parameters/observations
typedef coretools::Positive<double> TypeCovariateEffort;
typedef coretools::WeakType<double> TypeCovariateDetection;
typedef coretools::UnsignedInt<uint64_t> TypeCounts;
typedef coretools::UnsignedInt<uint64_t> TypeNu;
typedef coretools::WeakType<double> TypeGamma;
typedef coretools::Positive<double> TypeEffort;
typedef coretools::UnsignedInt<uint16_t> TypeTime;
typedef coretools::Probability TypeMuOrN;
typedef coretools::StrictlyPositive<double> TypeB;
typedef coretools::StrictlyPositiveWithMax<double, 2> TypeAlpha;
typedef coretools::WeakType<double> TypeBeta;
typedef coretools::WeakTypeWithExp<double> TypeLogPhi;
typedef coretools::StrictlyPositive<double> TypePhi;
typedef coretools::MinMaxVariable<double, 2> TypeLogSigma;

#endif // BIRPTYPES_H
