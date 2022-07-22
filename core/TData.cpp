//
// Created by caduffm on 7/20/22.
//

#include "TData.h"

using namespace coretools::instances;

bool TData::_shouldResizeEffort(size_t NumTimePointsCounts, size_t NumTimePointsEfforts) {
	return NumTimePointsCounts > 1 && NumTimePointsEfforts == 1;
}

void TData::_resizeEffort(size_t NumTimePointsCounts, coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts,
                          const std::shared_ptr<coretools::TNamesEmpty> &NamesTimePoints) {
	// expand effort to match dimensions of counts
	const size_t numMethods = Efforts.dimensions()[0];
	coretools::TMultiDimensionalStorage<TypeEffort, 2> newEffort;
	newEffort.prepareFillData(numMethods, {NumTimePointsCounts});

	for (size_t i = 0; i < numMethods; i++) {
		for (size_t j = 0; j < NumTimePointsCounts; j++) { newEffort.emplace_back(Efforts(i, 0)); }
	}

	newEffort.finalizeFillData();
	auto namesMethods = Efforts.getDimensionName(0);
	Efforts           = newEffort;

	// set dimensions names
	Efforts.setDimensionName(namesMethods, 0);
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
			_resizeEffort(NumTimePointsCounts, Efforts, NamesTimePoints);
			coretools::instances::logfile().list("Assuming constant efforts across all time points.");
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

void TData::_filterMethods(coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                           coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	// prepare new storage
	coretools::TMultiDimensionalStorage<TypeCounts, 2> filteredCounts;
	coretools::TMultiDimensionalStorage<TypeEffort, 2> filteredEfforts;
	filteredCounts.prepareFillData(Counts.dimensions()[0], {Counts.dimensions()[1]});
	filteredEfforts.prepareFillData(Efforts.dimensions()[0], {Efforts.dimensions()[1]});
	auto filteredMethodNames = std::make_shared<coretools::TNamesStrings>();

	for (size_t i = 0; i < Counts.dimensions()[0]; i++) {
		size_t hasEffort = 0;
		for (size_t k = 0; k < Counts.dimensions()[1]; k++) {
			if (Efforts(i, k) > 0.0) {
				hasEffort++;
			} else if (Efforts(i, k) == 0.0 && Counts(i, k) > 0) {
				UERROR("Invalid data for method ", (*Counts.getDimensionName(0))[i], " and time point ",
				       (*Counts.getDimensionName(1))[k], ": Counts can not be > 0 if effort is zero!");
			}
		}
		if (hasEffort >= 2) {
			// keep method counts and efforts
			filteredMethodNames->addName({(*Counts.getDimensionName(0))[i]});
			for (size_t k = 0; k < Counts.dimensions()[1]; k++) {
				filteredCounts.emplace_back(Counts(i, k));
				filteredEfforts.emplace_back(Efforts(i, k));
			}
		} else {
			logfile().warning("Removing method ", (*Counts.getDimensionName(0))[i],
			                  " because less than two time points have effort > 0!");
		}
	}

	filteredCounts.finalizeFillData();
	filteredEfforts.finalizeFillData();

	// check if any data are left
	if (filteredEfforts.dimensions()[0] == 0) { UERROR("No data set with sufficient effort!"); }

	// set dimension names
	filteredCounts.setDimensionName(filteredMethodNames, 0);
	filteredEfforts.setDimensionName(filteredMethodNames, 0);

	filteredCounts.setDimensionName(Counts.getDimensionName(1), 1);
	filteredEfforts.setDimensionName(Efforts.getDimensionName(1), 1);

	// overwrite storage
	Counts  = filteredCounts;
	Efforts = filteredEfforts;
}

void TData::read(coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                 coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	// read counts
	std::string filenameCounts = parameters().getParameterFilename("counts");
	_read(Counts, filenameCounts, "observational counts");
	const size_t numMethodsCounts    = Counts.dimensions()[0];
	const size_t numTimePointsCounts = Counts.dimensions()[1];

	// read efforts
	std::string filenameEff = parameters().getParameterFilename("efforts");
	_read(Efforts, filenameEff, "efforts");
	const size_t numMethodsEfforts    = Efforts.dimensions()[0];
	const size_t numTimePointsEfforts = Efforts.dimensions()[1];

	// check if dimensions match
	_checkDimensions_Methods(numMethodsCounts, numMethodsEfforts);
	_checkDimensions_TimePoints(numTimePointsCounts, numTimePointsEfforts, Efforts, Counts.getDimensionName(1));

	// check if dimension names match (i.e. same names for methods and time points)
	_checkDimensionNames(Counts, Efforts);

	// finally remove methods where < 2 time points have data
	_filterMethods(Counts, Efforts);
}

void TData::write(const std::string &Prefix, coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
                  coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts) {
	_write(Counts, Prefix + "_counts.txt", "observational counts");
	_write(Efforts, Prefix + "_efforts.txt", "efforts");
}