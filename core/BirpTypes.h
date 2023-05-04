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
typedef coretools::UnsignedInt<uint16_t> TypeIndex;
typedef coretools::WeakType<double> TypeGamma;
typedef coretools::Positive<double> TypeEffort;
typedef coretools::UnsignedInt<uint16_t> TypeTime;
typedef coretools::StrictlyPositiveWithMax<double> TypeMu;
typedef coretools::StrictlyPositive<double> TypeSigmaMu;
typedef coretools::MinMaxVariable<double> TypeLogA;
typedef coretools::StrictlyPositiveWithMax<double> TypeAlpha;
typedef coretools::WeakType<double> TypeBeta;
typedef coretools::StrictlyPositive<double> TypeSigmaBeta;
typedef coretools::WeakTypeWithExp<double> TypeLogPhi;
typedef coretools::StrictlyPositive<double> TypePhi;
typedef coretools::StrictlyPositive<double> TypeSigma;

#endif // BIRPTYPES_H
