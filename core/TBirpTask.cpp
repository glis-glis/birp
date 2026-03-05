#include "TBirpCore.h"
#include "TBirpTask.h"

TTask_infer::TTask_infer() : coretools::TTask("Inferring rates of change in species densities from count data.") {}

void TTask_infer::run() {
	TBirpCore birp;
	birp.infer();
}

TTask_simulate::TTask_simulate() : coretools::TTask("Simulating count data under Birp model.") {}

void TTask_simulate::run() {
	TBirpCore birp;
	birp.simulate();
}
