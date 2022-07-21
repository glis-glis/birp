//
// Created by caduffm on 7/20/22.
//

#include "TData.h"

using namespace coretools::instances;

bool TData::_shouldResizeEffort(size_t NumTimePointsCounts, size_t NumTimePointsEfforts) {
	return NumTimePointsCounts > 1 && NumTimePointsEfforts == 1;
}

void TData::_resizeEffort(size_t NumTimePointsCounts, size_t NumMethods,
                          coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts,
                          const std::shared_ptr<coretools::TNamesEmpty> &NamesTimePoints) {
	// expand effort to match dimensions of counts
	coretools::TMultiDimensionalStorage<TypeEffort, 2> newEffort;
	newEffort.prepareFillData(NumMethods, {NumTimePointsCounts});

	for (size_t i = 0; i < NumMethods; i++) {
		for (size_t j = 0; j < NumTimePointsCounts; j++) { newEffort.emplace_back(Efforts(i, 0)); }
	}

	newEffort.finalizeFillData();
	Efforts = newEffort;

	// set time point names from data
	Efforts.setDimensionName(NamesTimePoints, 1);
}

void TData::_checkDimensions_Methods(size_t NumMethodsCounts, size_t NumMethodsEfforts) {
	// check if dimensions match: number of methods
	if (NumMethodsCounts != NumMethodsEfforts) {
		UERROR("Detected a different number of methods in counts file (", NumMethodsCounts, ") than in efforts file (",
		       NumMethodsEfforts, ")!");
	}
}

void TData::_checkDimensions_TimePoints(size_t NumTimePointsCounts, size_t NumTimePointsEfforts,
                                        coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts,
                                        const std::shared_ptr<coretools::TNamesEmpty> &NamesTimePoints) {
	// check if dimensions match: number of time points
	if (NumTimePointsCounts != NumTimePointsEfforts) {
		if (_shouldResizeEffort(NumTimePointsCounts, NumTimePointsEfforts)) {
			// Counts were measured at >1 time points, but effort is constant across all time points -> resize effort
			_resizeEffort(NumTimePointsCounts, NumTimePointsEfforts, Efforts, NamesTimePoints);
		} else {
			UERROR("Detected a different number of time points in counts file (", NumTimePointsCounts,
			       ") than in efforts file (", NumTimePointsEfforts, ")!");
		}
	}
}

void TData::_checkDimensionNames(const coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                                 const coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	if (*Counts.getDimensionName(0) != *Efforts.getDimensionName(0)) {
		UERROR("Methods (row names) in the counts file do not match the methods in the efforts file!");
	}
	if (*Counts.getDimensionName(1) != *Efforts.getDimensionName(1)) {
		UERROR("Time points (column names) in the counts file do not match the time points in the efforts file!");
	}
}

void TData::read(coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                 coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	// read counts
	std::string filenameCounts = parameters().getParameterFilename("counts");
	_read(Counts, filenameCounts, "observational counts");
	const size_t numMethodsCounts    = Counts.dimensions()[0];
	const size_t numTimePointsCounts = Counts.dimensions()[1];

	// read efforts
	std::string filenameEff = parameters().getParameterFilename("effort");
	_read(Efforts, filenameEff, "efforts");
	const size_t numMethodsEfforts    = Efforts.dimensions()[0];
	const size_t numTimePointsEfforts = Efforts.dimensions()[1];

	// check if dimensions match
	_checkDimensions_Methods(numMethodsCounts, numMethodsEfforts);
	_checkDimensions_TimePoints(numTimePointsCounts, numTimePointsEfforts, Efforts, Counts.getDimensionName(1));

	// check if dimension names match (i.e. same names for methods and time points)
	_checkDimensionNames(Counts, Efforts);
}

void TData::write(const std::string &Prefix, coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                  coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	_write(Counts, Prefix + "_counts.txt", "observational counts");
	_write(Efforts, Prefix + "_efforts.txt", "efforts");
}