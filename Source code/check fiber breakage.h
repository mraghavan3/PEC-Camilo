#include "classes and structures.h"
#include "fiber.h"
#include "Broadphase.h"

void checkFibersBreakage(vector<fiber> & fibers, double criticalAngle,thrust::host_vector<int> &fiberIndices, thrust::host_vector<int> &segmentIndices );