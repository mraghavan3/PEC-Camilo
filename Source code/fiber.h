#pragma once
#include "classes and structures.h"
#include "segment.h"
class fiber
{
public:
	
	int numberOfHinges;
	int numberOfSegments;
	int firstHinge;
	vector<Hinge> hinges;
	vector<Vector3d> midpoints;
	vector<double> segmentsLength;
	bool isBroken;
	int cellID;
	vector<segment> segments;
	int elementID;
	fiber splitFiber(int hingeID);

};
