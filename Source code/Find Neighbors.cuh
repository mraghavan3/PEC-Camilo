#pragma once
#include "classes and structures.h"
#include "fiber.h"
#include "utilities.h"
#include "Broadphase.h"
struct Cell
{
	// each cell could contain a segment that can be accesed as a double index ( fiber index , segment index)

	vector<int>fiberIndex;
	vector<int>segmentIndex;
	vector<int>neighborFiberIndex;
	vector<int>neighborSegmentIndex;
	bool hasFibers;
	int a;
	
};


void findNeighbours( vector<fiber> &fibers, double fiberRadius, thrust::host_vector<long long> &possibleCollisions, thrust::host_vector<int> fiberIndices);
void findBoundingBoxes(vector<fiber> fibers, double * & boundingBoxes, double fiberRadius); 