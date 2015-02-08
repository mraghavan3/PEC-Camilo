#include "Find Neighbors.cuh"
#include "Broadphase.h"

void getMaxMin( double a, double b , double &max , double &min ){

	if (a > b )
			{
				max = a;
				min = b;
			}
			else
			{
				min = a;
				max = b;
			}

}

void appendVector( vector<int> & destinationVector, vector<int> sourceVector){
	int initialSize = destinationVector.size();
	destinationVector.resize(initialSize + sourceVector.size()); // added to improve efficiency
	for (uint i = 0; i < sourceVector.size(); i++)
	{
		destinationVector[initialSize + i ] = sourceVector[i];
		//destinationVector.push_back(sourceVector[i]); //removed to improve efficiency
	}
}


void findNeighbours( vector<fiber> &fibers, double fiberRadius, thrust::host_vector<long long> & potentialCollisions, thrust::host_vector<int> fiberIndices){
	vector<Cell> cells ;
	
	
	int timeA = getMilliCount();

	// first step, find box dimensions, and cell dimension rc
	// I don't like having to do this, it seems like alot of operations to find rc
	
	double maxSegmentLength=0;
	double minX  = fibers[0].hinges[0].position(0);
	double minY = fibers[0].hinges[0].position(1);
	double minZ = fibers[0].hinges[0].position(2);
	double maxX = fibers[0].hinges[0].position(0);
	double maxY = fibers[0].hinges[0].position(1);
	double maxZ = fibers[0].hinges[0].position(2);
	
	for (uint i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j <fibers[i].numberOfSegments ; j++)
		{
			//compute length of each segment
			
			double segmentLength = sqrt(  (fibers[i].hinges[j].position-fibers[i].hinges[j+1].position).dot(fibers[i].hinges[j].position-fibers[i].hinges[j+1].position));
			
			// record the longest segment
			if (segmentLength>maxSegmentLength){maxSegmentLength = segmentLength;}
			
			// compute midpoint of current segment
			
			fibers[i].midpoints[j] = (fibers[i].hinges[j].position + fibers[i].hinges[j+1].position)/2;
			
			//record the minimum and maximum coords of the fibers
			if (minX > fibers[i].midpoints[j](0)){minX = fibers[i].midpoints[j](0);}
			if (minY > fibers[i].midpoints[j](1)){minY = fibers[i].midpoints[j](1);}
			if (minZ > fibers[i].midpoints[j](2)){minZ = fibers[i].midpoints[j](2);}
			if (maxX < fibers[i].midpoints[j](0)){maxX = fibers[i].midpoints[j](0);}
			if (maxY < fibers[i].midpoints[j](1)){maxY = fibers[i].midpoints[j](1);}
			if (maxZ < fibers[i].midpoints[j](2)){maxZ = fibers[i].midpoints[j](2);}

		}
	}
	
	// this is the characteristic length of the cell 
	double rc =maxSegmentLength;
	
	// these are the dimensions of the box
	Vector3d boxLength;//L
	boxLength<< (maxX-minX), (maxY - minY) , (maxZ-minZ) ;
	
	//cout<<endl<< "box dimensions" << boxLength;
	//cout<<endl<< "characteristic dimension " << rc;
	//cout<<endl<< " minx " << minX << " maxx " << maxX << endl;
	
	// the number of cells in each direction
	Vector3i cellsNumber; // Lc
	cellsNumber<< (int)(floor(boxLength(0)/rc)) , (int)(floor(boxLength(1)/rc)) , (int)(floor(boxLength(2)/rc)) ;
	for (int i = 0; i < 3; i++)
	{
		if (cellsNumber(i) ==0)
		{
			cellsNumber(i) =1;
		}
	}


	// 
	Vector3d cellLength; // rc
	cellLength<< (boxLength(0)/cellsNumber(0)) , (boxLength(1)/cellsNumber(1)) , (boxLength(2)/cellsNumber(2)) ;
	/*

	// 
	int totalNumberCells = cellsNumber(0)*cellsNumber(1)*cellsNumber(2);
	//cout<<" totalNumber of cells" << totalNumberCells<<endl;
	cells.resize(totalNumberCells); // For a reason I cant understand I have to have this line and the for that follows, 
									// otherwise I'll get an error that I do not know what does it mean;
	//cout<< endl << "number of cells " << cellsNumber(0)<< " y " <<cellsNumber(1) << " z " <<cellsNumber(2)<<endl;
	
	//cout<<"Number of cells " << cells.size() << " Number of fibers " << fibers.size()<< endl;


	for (int i = 0; i < totalNumberCells+1; i++)
	{
		Cell tempcell;
		tempcell.hasFibers =false;
		//cells[i] = tempcell;
		tempcell.fiberIndex.clear();

		cells.push_back(tempcell);
		cells[i].hasFibers=false;
		
	}

	 // Lets find the cell where each segment is 

	for (uint i = 0; i < fibers.size(); i++)
	{
		//cout<< "fiber " << i << " of " << fibers.size() << endl;
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			
			Vector3i cellID;
			cellID(0) = floor( (fibers[i].midpoints[j].x() - minX) / cellLength.x()); 
			cellID(1) = floor( (fibers[i].midpoints[j].y() - minY) / cellLength.y()); 
			cellID(2) = floor( (fibers[i].midpoints[j].z() - minZ) / cellLength.z()); 

			//cout<< " index of segment x " << cellID.x() << " y " <<cellID.y() << " z "<< cellID.z() << endl;
			int cellIndex = cellID.x() * cellsNumber.y() * cellsNumber.z() + cellID.y() * cellsNumber.z() + cellID.z();
			//cout<< "segment is on cell number : " << cellIndex << " max index "<<  totalNumberCells  << endl; 
			cells[cellIndex].hasFibers =true;
			//cout<< " done with this cycle "<<endl;
			cells[cellIndex].fiberIndex.push_back(i); // I need to pushback because  I dont know the size
			//cells[cellIndex].a =2;
			cells[cellIndex].segmentIndex.push_back(j);
			fibers[i].segments[j].cellID = cellIndex;


			//
		}

	}
	

	// Lets find the neighbors of every cell
	// iterating over the three dimensional indices of the cells
	for (int xIndex = 0; xIndex < cellsNumber.x(); xIndex++)
	{
		for (int yIndex = 0; yIndex < cellsNumber.y(); yIndex++)
		{
			for (int zIndex = 0; zIndex < cellsNumber.z(); zIndex++)
			{
				// name of variable says it all
				int currentCellIndex = xIndex * cellsNumber.y() * cellsNumber.z() + yIndex * cellsNumber.z() + zIndex;
				//cout<< " current index " << currentCellIndex<<endl;
				int numberOFneighboringCells =0 ;
				if ( cells[currentCellIndex].hasFibers){
					

				// Now I have to look for the 26 neighboring cells

					for (int xNeighbor = xIndex-1; xNeighbor < xIndex+2; xNeighbor++)
					{
						if (xNeighbor >= 0 && xNeighbor < cellsNumber.x()) // the index needs to be bigger than or equal to zero
																		   // And smaller than the number of cell in x direction
						{
							for (int yNeighbor = yIndex-1; yNeighbor < yIndex+2; yNeighbor++)
							{
								if (yNeighbor >= 0 && yNeighbor < cellsNumber.y())
								{
									for (int zNeighbor = zIndex-1; zNeighbor < zIndex + 2; zNeighbor++)
									{
										if (zNeighbor >= 0 && zNeighbor < cellsNumber.z())
										{
											int neighborCellIndex = xNeighbor * cellsNumber.y() * cellsNumber.z() + yNeighbor * cellsNumber.z() + zNeighbor;
											if (cells[neighborCellIndex].hasFibers)
											{
												//copy fiber index vector and segment index vector of the neighbor cell to the neighbor fiber array
												//and neighbor index vector
												appendVector(cells[currentCellIndex].neighborFiberIndex , cells[neighborCellIndex].fiberIndex);
												appendVector(cells[currentCellIndex].neighborSegmentIndex , cells[neighborCellIndex].segmentIndex);
											}
											
											if (neighborCellIndex!= currentCellIndex)
											{
												numberOFneighboringCells++;
											}
										}
									}
								}
							}
						}
					}
				}//ends first if
				//cout<< " number of neighboring cells " << numberOFneighboringCells << endl;
			}
		}
	}

	//cout<< " total number of cells " << totalNumberCells; 

	
 // copy the list of neighbors to each segment

	for (uint i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
		 fibers[i].segments[j].neighboringFibers.resize(0);
		 fibers[i].segments[j].neighboringSegment.resize(0);
		 appendVector(fibers[i].segments[j].neighboringFibers, cells[fibers[i].segments[j].cellID].neighborFiberIndex);
		 appendVector(fibers[i].segments[j].neighboringSegment, cells[fibers[i].segments[j].cellID].neighborSegmentIndex);
		}
	}

	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
	cout << "time required to find neighbors my method " << getMilliSpan(timeA)<<" ms"<<endl;
	#endif


	*/

	// Broad phase use
	
	
	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
				 timeA =getMilliCount();
	#endif

	// Step 1, get AABB data

	thrust::host_vector<real3> aabb_data_H;
	//custom_vector<real3> aabb_data;
	potentialCollisions.clear();

	double offsetFactor=1.5;


	// first half of the array has the minimum point of the bounding box, that part will be filled with the followin nested loops
	for (uint i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			double minCoord;
			double maxCoord;
			real3 temp;

			getMaxMin(fibers[i].hinges[j].position.x(),fibers[i].hinges[j+1].position.x() ,maxCoord,minCoord);
			temp.x = minCoord - offsetFactor*fiberRadius;
			getMaxMin(fibers[i].hinges[j].position.y(),fibers[i].hinges[j+1].position.y() ,maxCoord,minCoord);
			temp.y = minCoord - offsetFactor*fiberRadius;
			getMaxMin(fibers[i].hinges[j].position.z(),fibers[i].hinges[j+1].position.z() ,maxCoord,minCoord);
			temp.z = minCoord - offsetFactor*fiberRadius;
			aabb_data_H.push_back(temp);
			

		}
	}
	for (uint i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			double minCoord;
			double maxCoord;
			real3 temp;

			getMaxMin(fibers[i].hinges[j].position.x(),fibers[i].hinges[j+1].position.x() ,maxCoord,minCoord);
			temp.x = maxCoord + offsetFactor*fiberRadius;
			getMaxMin(fibers[i].hinges[j].position.y(),fibers[i].hinges[j+1].position.y() ,maxCoord,minCoord);
			temp.y = maxCoord + offsetFactor*fiberRadius;
			getMaxMin(fibers[i].hinges[j].position.z(),fibers[i].hinges[j+1].position.z() ,maxCoord,minCoord);
			temp.z = maxCoord + offsetFactor*fiberRadius;
			aabb_data_H.push_back(temp);

		}
	}

	#ifdef SIM_ENABLE_GPU_MODE
		custom_vector<real3> aabb_data = aabb_data_H;
	#else
		
	#endif
	
	
	//Step 2 run broadphase algorithm
	
	Broadphase broadphaseManager;
	broadphaseManager.setBinsPerAxis(make_real3(cellsNumber.x(),cellsNumber.y(),cellsNumber.z()));

	#ifdef SIM_ENABLE_GPU_MODE
		custom_vector<long long> potentialCollisions_D  = potentialCollisions;
		broadphaseManager.detectPossibleCollisions(aabb_data, potentialCollisions_D);
		potentialCollisions = potentialCollisions_D;
	#else
		broadphaseManager.detectPossibleCollisions(aabb_data_H, potentialCollisions);
	#endif
	


	/*cout<<"Number of possible contacts"<< broadphaseManager.getNumPossibleContacts()<<endl;
	cout<< " number of neighboring fibers" << fibers[0].segments[0].neighboringFibers.size()<< endl;
	cout<< " number of neighboring fibers" << fibers[0].segments[0].neighboringSegment.size()<< endl;
	//cout << "time required to find neighbors with broadphase 1 " << getMilliSpan(timeA)<<" ms"<<endl;
	//Remove collisions beteen segments of the same fiber*/

	thrust::host_vector<int> stencil_H;
	
	
	for (uint i = 0; i < broadphaseManager.getNumPossibleContacts(); i++)
	{
		long particleA = ((long)(potentialCollisions[i]>>32));
		long particleB = (long) potentialCollisions[i];
		if (fiberIndices[particleA] == fiberIndices[particleB])
		{
			stencil_H.push_back(1);
			//cout<<" Same fiber neighbors " <<endl;
		}
		else{
			stencil_H.push_back(0);
		}

	}

	custom_vector<int> stencil= stencil_H;

	#ifdef PRINT_DEBUG_INFO
	cout<<"number of possible collisions " << potentialCollisions.size()<<endl;
	cout << "time required to find neighbors with broadphase 1 " << getMilliSpan(timeA)<<" ms"<<endl;
	#endif
/*	for (int i = 0; i < broadphaseManager.getNumPossibleContacts(); i++)
	{


		//cout<<"My Method "<<" fiber "<<fibers[0].segments[2].neighboringFibers[i] <<" Segment "<< fibers[0].segments[2].neighboringSegment[i]<<endl;

		int fiberIndexA = fiberIndices[((long)(potentialCollisions[i]>>32)) ];
		int segmentA = ((long)(potentialCollisions[i]>>32)) - segmentIndex[fiberIndexA];
		int fiberIndexB = fiberIndices[(long) potentialCollisions[i]];
		int segmentB = (long) potentialCollisions[i] - segmentIndex[fiberIndexB];

		cout<<"BroadPhase,  element A  " <<  ((long)(potentialCollisions[i]>>32)) << " element B " << (long) potentialCollisions[i]<< " stencil: " << stencil[i]<<endl;
		cout<<" fiber A " << fiberIndexA << " segment " << segmentA <<" Fiber B " <<fiberIndexB;
				cout<<" segment " << segmentB<<endl;
	}*/


	int numberOfCollisions;

	#ifdef SIM_ENABLE_GPU_MODE
	//custom_vector<long long> potentialCollisions_D  = potentialCollisions;
	if (potentialCollisions_D.size() != 0)
	{
		numberOfCollisions = broadphaseManager.removeSameFiberCollisions(potentialCollisions_D,stencil);
		potentialCollisions = potentialCollisions_D;

		int realCollissions = stencil_H.size() - numberOfCollisions;
		potentialCollisions.resize(realCollissions);

		/*
		cout<<"Stencil : ";
		for (int i = 0; i < stencil_H.size(); i++)
		{
			cout<<stencil_H[i] << ", ";
		}
		cout<<endl;*/
	}
		
		
	#else
	if (potentialCollisions.size() != 0 )
	{
		numberOfCollisions = broadphaseManager.removeSameFiberCollisionsH(potentialCollisions,stencil);
	}
		
	#endif

	
	
	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
	 cout << "time required to find neighbors with broadphase" << getMilliSpan(timeA)<<" ms"<<endl;
	#endif



}
       





void findBoundingBoxes( vector<fiber> fibers, double * &boundingBox, double fiberRadius){

	// bounding box is an array with the format:
	// minx miny minz maxx maxy maxz

	double minx, miny , minz, maxx, maxy, maxz;
	int index =0;

	// loops over every hinge of the system
	for (uint i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			
			// compares each pair of coordinate components to find out whis is the minimum and maximum
			// 
			getMaxMin(fibers[i].hinges[j].position.x() , fibers[i].hinges[j+1].position.x(), maxx , minx); 
			getMaxMin(fibers[i].hinges[j].position.y() , fibers[i].hinges[j+1].position.y(), maxy , miny);
			getMaxMin(fibers[i].hinges[j].position.z() , fibers[i].hinges[j+1].position.z(), maxx , minz);


			// stores the data
			boundingBox[index] = minx - fiberRadius; index++;
			boundingBox[index] = miny - fiberRadius; index++;
			boundingBox[index] = minz - fiberRadius; index++;
			boundingBox[index] = maxx + fiberRadius; index++;
			boundingBox[index] = maxy + fiberRadius; index++;
			boundingBox[index] = maxz + fiberRadius; index++;
		}
	}


}
