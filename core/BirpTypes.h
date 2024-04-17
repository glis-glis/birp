//
// Created by caduffm on 7/20/22.
//

#ifndef BIRPTYPES_H
#define BIRPTYPES_H

#include "coretools/Types/commonWeakTypes.h"
#include "coretools/Types/probability.h"
#include <cstdint>

//-------------------------------
// type alias
//-------------------------------

// gamma
using TypeGamma                     = coretools::Unbounded;
constexpr static size_t NumDimGamma = 1;

// time of change
using TypeTime                             = coretools::Positive;
constexpr static size_t NumDimTimeOfChange = 1;

// alpha
using TypeAlpha                     = coretools::StrictlyPositiveWithMax<0>;
constexpr static size_t NumDimAlpha = 1;

// beta
using TypeBeta                     = coretools::Unbounded;
constexpr static size_t NumDimBeta = 1;

// muOrN
using TypeMuOrN                     = coretools::Probability;
constexpr static size_t NumDimMuOrN = 1;

// b
using TypeB                     = coretools::StrictlyPositive;
constexpr static size_t NumDimB = 1;

// logPhi
using TypeLogPhi                     = coretools::Unbounded;
constexpr static size_t NumDimLogPhi = 2;

// logSigma
using TypeLogSigma                     = coretools::MinMaxVariable<1>;
constexpr static size_t NumDimLogSigma = 1;

// counts
using TypeCounts                     = coretools::UnsignedInt64;
constexpr static size_t NumDimCounts = 2;

// other types
using TypeCovariateEffort    = coretools::Positive;
using TypeCovariateDetection = coretools::Unbounded;
using TypeNu                 = coretools::UnsignedInt64;
using TypeEffort             = coretools::Positive;
using TypePhi                = coretools::StrictlyPositive;

#endif // BIRPTYPES_H
