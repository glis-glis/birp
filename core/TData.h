//
// Created by caduffm on 7/20/22.
//

#ifndef TDATA_H
#define TDATA_H

#include "BirpTypes.h"
#include "TDataFile.h"
#include "TLog.h"
#include "TParameters.h"

class TData {
private:
public:
	static void read(coretools::TMultiDimensionalStorage<TypeCounts, 2> &Data);
	static void write(const std::string &Prefix, coretools::TMultiDimensionalStorage<TypeCounts, 2> &Data);
};

#endif // TDATA_H
