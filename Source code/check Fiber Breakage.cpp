#include "check fiber breakage.h"


double dotProduct_(double Ax, double Ay, double Az, double Bx, double By, double Bz  ){

	return Ax*Bx + Ay*By + Az*Bz;

}

void checkFibersBreakage(vector<fiber> &fibers,double criticalAngle,thrust::host_vector<int> &fiberIndices, thrust::host_vector<int> &segmentIndices){
	bool fiberBroke = false;
	int numInFibers= fibers.size();
	for (uint i = 0; i < numInFibers; i++) //for each fiber
	{
		for (int j = 1; j < fibers[i].numberOfSegments; j++) // for each hinge in each fiber
		{
			double angle = acos((fibers[i].hinges[j-1].position - fibers[i].hinges[j].position).dot(fibers[i].hinges[j+1].position - fibers[i].hinges[j].position )    );
			//cout<<" angle " <<angle<< " critical angle  " <<criticalAngle << " " ;
			if (angle <= criticalAngle)
			{


				cout<< " fiber breaks";
				fibers.push_back(fibers[i].splitFiber(j));
				fiberBroke = true;
			}

		}

	}

	if (fiberBroke)
	{
		segmentIndices.clear();
		fiberIndices.clear();
		// We need to regenerate the fiber indices and segment indices.
		int count=0;
	
		for (uint i = 0; i < fibers.size(); i++)
		{
			segmentIndices.push_back(count);
			for (int j = 0; j < fibers[i].numberOfSegments; j++)
			{
				fiberIndices.push_back(i);
				count++;
			}
		
		}

	}


}