#ifndef TBIRPTASK_H_
#define TBIRPTASK_H_

#include "coretools/Main/TTask.h"

class TTask_infer : public coretools::TTask {
public:
	TTask_infer();
	void run() override;
};

class TTask_simulate : public coretools::TTask {
public:
	TTask_simulate();
	void run() override;
};

#endif
