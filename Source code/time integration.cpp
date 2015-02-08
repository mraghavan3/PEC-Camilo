#include "time integration.h"

int updatePositions(vector<fiber> &fibers, double dt){

	// Explicit Euler Scheme

	// we want to preserve the length of the segments, that is way the lengths are used

	for (int i = 0; i < fibers.size(); i++) // loop through all the fibers
	{
		// if the fiber is stationary dont do anything
		if (fibers[i].hinges[0].isStationary)
		{

		}else
		{



			// Stores the length of the first segment, when the second hinge is updated
			// the distance between the first and second will be the same
			Vector3d vectorialDistanceOld = fibers[i].hinges[1].position -fibers[i].hinges[0].position;
			double oldLength = sqrt(vectorialDistanceOld.dot(vectorialDistanceOld));


			fibers[i].hinges[0].position +=  fibers[i].hinges[0].velocity*dt; // update first hinge of the fiber

			for (int j = 1; j < fibers[i].numberOfHinges; j++) // loop through all the hinges of each fiber
			{

				Vector3d vectorialDistanceNextHinge ;
				double segmentLength =0  ;

				// stores length of the next segment
				if ( j < fibers[i].numberOfHinges-1)
				{
					vectorialDistanceNextHinge = fibers[i].hinges[j+1].position -fibers[i].hinges[j].position;
					segmentLength = sqrt(vectorialDistanceNextHinge.dot(vectorialDistanceNextHinge));
				}

				// update the current hinge

				fibers[i].hinges[j].position += fibers[i].hinges[j].velocity *dt;

				// find the orientation of the segment
				Vector3d segmentOrientation = fibers[i].hinges[j].position - fibers[i].hinges[j-1].position;
				segmentOrientation = segmentOrientation/(sqrt(segmentOrientation.dot(segmentOrientation)) ); // Normalized orientation 

				// adjust the length of the segment
				fibers[i].hinges[j].position = fibers[i].hinges[j-1].position + segmentOrientation * oldLength;

				oldLength = segmentLength;

			}
		}
	}

	

	return 0;
}