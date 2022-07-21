//
// Created by caduffm on 7/20/22.
//

#ifndef TDATA_H
#define TDATA_H

#include "BirpTypes.h"
#include "TDataFile.h"
#include "TError.h"
#include "TLog.h"
#include "TNames.h"
#include "TParameters.h"

// TODO: Remove methods with less than two time points with data (effort > 0)

class TData {
private:
	template<typename T>
	static void _read(coretools::TMultiDimensionalStorage<T, 2> &Storage, const std::string &Filename,
	                  const std::string &Description) {
		// both counts and efforts have the same file format, but different types -> templated it
		using namespace coretools::instances;

		// prepare objects for storing row- and column names
		auto methodNames = std::make_shared<coretools::TNamesStrings>();
		auto timeNames   = std::make_shared<coretools::TNamesStrings>();

		// guess memory: very few methods, doesn't matter
		const size_t guessNumMethods = 2;

		// read!
		logfile().listFlush("Reading ", Description, " from file '" + Filename + "'...");
		coretools::TDataReaderBase<T, 2> reader;
		reader.read(Filename, &Storage, {methodNames, timeNames}, {'\n', '\t'}, {true, true}, {0, 0}, {false, false},
		            guessNumMethods, coretools::colNames_multiline);
		logfile().done();

		// report
		logfile().conclude("Detected ", Storage.dimensions()[0], " methods and ", Storage.dimensions()[1],
		                   " time points.");
	}

	template<typename T>
	static void _write(coretools::TMultiDimensionalStorage<T, 2> &Storage, const std::string &Filename,
	                   const std::string &Description) {
		// both counts and efforts have the same file format, but different types -> templated it
		using namespace coretools::instances;

		logfile().listFlush("Writing ", Description, " to file '" + Filename + "'...");
		coretools::TDataWriterBase<T, 2> writer;
		writer.write(Filename, &Storage, {'\n', '\t'}, {true, true}, coretools::colNames_multiline);
		logfile().done();
	}

	static bool _shouldResizeEffort(size_t NumTimePointsCounts, size_t NumTimePointsEfforts);
	static void _resizeEffort(size_t NumTimePointsCounts, coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts,
	                          const std::shared_ptr<coretools::TNamesEmpty> &NamesTimePoints);
	static void _checkDimensions_Methods(size_t NumMethodsCounts, size_t NumMethodsEfforts);
	static void _checkDimensions_TimePoints(size_t NumTimePointsCounts, size_t NumTimePointsEfforts,
	                                        coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts,
	                                        const std::shared_ptr<coretools::TNamesEmpty> &NamesTimePoints);
	static void _checkDimensionNames(const coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
	                                 const coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts);

public:
	static void read(coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
	                 coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts);
	static void write(const std::string &Prefix, coretools::TMultiDimensionalStorage<TypeCounts, 2> &Counts,
	                  coretools::TMultiDimensionalStorage<TypeEffort, 2> &Efforts);
};

#endif // TDATA_H
