#include "Excluded Volume.h"


struct Distance
{
	Vector3d Gab;
	Vector3d pb;
	double Gab_norm;
	Vector3d r;
	Vector3d r2;
	bool collisionCourse;

};




void Distance_point_point(Vector3d pointA, Vector3d pointB, Vector3d & g_ab, double & g_ab_norm){

	g_ab = pointA-pointB;
	g_ab_norm = sqrt(g_ab.dot(g_ab));
	
}

void Distance_segments(Vector3d hingeA0, Vector3d hingeA1, Vector3d hingeB0, Vector3d hingeB1, Vector3d & r, Vector3d & r2,Vector3d & g_ab, double & g_ab_norm, bool & collisionCourse  ){
	/*Vector3d pa = hingeA1 - hingeA0;
	Vector3d pb = hingeB1 - hingeB0; 
	

	double la = sqrt(pa.dot(pa));
	double lb = sqrt(pb.dot(pb));


	pa = pa/la;
	pb = pb/lb;
	double padotpa = pa.dot(pa);
	double pbdotpb = pb.dot(pb);
	double padotpb = pa.dot(pb);
	Vector3d vecDifference = hingeA0-hingeB0;

	
	double Sab = ( (vecDifference).dot(pb) * (padotpb- (vecDifference).dot(pa) ) ) / ( 1 - padotpb* padotpb ) ;
	double Sba = ( (-vecDifference).dot(pa) * (padotpb - (-vecDifference).dot(pb) ) ) / ( 1 - padotpb* padotpb ) ;
	

	/*
	double Sab = ( (hingeA0-hingeB0).dot(pb) * ( pa.dot(pb) - (hingeA0-hingeB0).dot(pa) ) ) / ( 1 - pa.dot(pb)* pa.dot(pb) ) ;
	double Sba = ( (hingeB0-hingeA0).dot(pa) * ( pb.dot(pa) - (hingeB0-hingeA0).dot(pa) ) ) / ( 1 - pb.dot(pa)* pb.dot(pa) ) ;
	
	r= pb * Sba;
	r2= pa * Sba;

	g_ab = hingeA0 + Sab*pa -hingeB0 - Sba*pb;
	g_ab_norm = sqrt(g_ab.dot(g_ab));
	collisionCourse =false;

	if (Sab <= la && Sab >= 0){
		if ( Sba <= lb && Sba >= 0){
			collisionCourse =true;
		}
	}
	*/
	// New test
	Vector3d u = hingeA1 - hingeA0;
	Vector3d v = hingeB1 - hingeB0;
	Vector3d w = hingeA0 - hingeB0;
	double a = u.dot(u);      // square of the length of segment A    
	double b = u.dot(v);
	double c = v.dot(v);      // square root of the length of segment B  
	double d = u.dot(w);
	double e = v.dot(w);
	double D = a*c - b*b;      

	double SMALL_NUM = 0.0000000001;
	double sc,tc;

    if (D < SMALL_NUM)  
	{
        sc = 0.0;
		if(b>c) {tc = d/b;} else {tc = e/c;} 
	}   
    else 
	{
        sc = (b*e - c*d) / D;
        tc = (a*e - b*d) / D;
	} 
    g_ab =( w + (sc * u) - (tc * v));
    g_ab_norm = sqrt(g_ab.dot(g_ab));

	r = -u * sc;
	r2 = -v * tc;
	if (sc <= 1 && sc >= 0){
		if ( tc <= 1 && tc >= 0){
			collisionCourse =true;
		}
	}
}

void Distance_point_segment(Vector3d rA , Vector3d rB, Vector3d  rB_end, Vector3d & r, Vector3d & g_ab, double & g_ab_norm, bool & collisionCourse ){

	Vector3d pb = rB_end - rB;
	double lb = sqrt(pb.dot(pb));
	pb = pb/lb;
	double Sba = (rA - rB).dot(pb);
	g_ab = rA -rB - Sba*pb;
	g_ab_norm = sqrt(g_ab.dot(g_ab));

	r= pb * Sba;
	collisionCourse =false;

	
		if ( Sba <= lb && Sba >= 0){
			collisionCourse =true;
		}
	

}




// compute excluded volume forces for debugging
void computeExludedVolumeForces( Hinge& HingeA0,  const Hinge& HingeA1, Hinge& HingeB0, const Hinge& HingeB1, double threshold, double Excl_vol_fac, double bead_radius, double & gab_min, vector<Distance> & Vec, bool isForNeighborsList, int timeStep){
	// TO DO, this is too verbose, make it more compact

	//vector<Distance> Vec(7);
	Distance tempVec;

	Distance_segments(HingeA0.position,HingeA1.position,HingeB0.position,HingeB1.position, tempVec.r, tempVec.r2, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	Vec[0]=tempVec; // 1

	Distance_point_segment(HingeA0.position, HingeB0.position, HingeB1.position, tempVec.r, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	tempVec.r2.setZero();
	Vec[1]=tempVec; // 2

	Distance_point_segment(HingeA1.position, HingeB0.position, HingeB1.position, tempVec.r, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	Vec[2]=tempVec; // 3

	Distance_point_point(HingeA0.position, HingeB0.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r.setZero();
	tempVec.r2.setZero();
	tempVec.collisionCourse =true; 
	Vec[3]=tempVec; // 4

	Distance_point_point(HingeA1.position, HingeB0.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r.setZero();
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	tempVec.collisionCourse =true;
	Vec[4]=tempVec; // 5

	Distance_point_point(HingeA0.position, HingeB1.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r = HingeB1.position-HingeB0.position;
	tempVec.r2.setZero();
	tempVec.collisionCourse =true;
	Vec[5]=tempVec; // 6

	Distance_point_point(HingeA1.position, HingeB1.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r = HingeB1.position-HingeB0.position;
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	tempVec.collisionCourse =true;
	Vec[6]=tempVec; // 7

	gab_min = Vec[6].Gab_norm;
	int k =6;

	for (int i = 0; i < 6; i++)
	{
		if (Vec[i].Gab_norm < gab_min  && Vec[i].collisionCourse)
		{
			gab_min = Vec[i].Gab_norm;
			k=i;

			


		}
	}
	
	if (!isForNeighborsList && gab_min < threshold && Vec[k].collisionCourse )
	{
		Vector3d Excluded_volume_partial_force = -Vec[k].Gab * Excl_vol_fac * exp(-2 *(Vec[k].Gab_norm / bead_radius -2));
		HingeB0.exluded_volume_force +=  Excluded_volume_partial_force ;
		HingeB0.excluded_volume_torque += Vec[k].r2.cross(Excluded_volume_partial_force);


		HingeA0.exluded_volume_force -=  Excluded_volume_partial_force ;
		HingeA0.excluded_volume_torque += Vec[k].r.cross(-Excluded_volume_partial_force);



		/*
		cout<<" segment A  : " << HingeA0.position.transpose() << " - " << HingeA1.position.transpose()<<endl;
		cout<<" segment B  : " << HingeB0.position.transpose() << " - " << HingeB1.position.transpose()<<endl;
		
		cout<<" segment A Excluded volume force : "<< HingeA0.exluded_volume_force.transpose() << endl;
		cout<<" segment B Excluded volume force : "<< HingeB0.exluded_volume_force.transpose() << endl;
		*/
	}

	if (timeStep% 1000 == 0)
	{
		cout<<" segment A  : " << HingeA0.position.transpose() << " - " << HingeA1.position.transpose()<<endl;
		cout<<" segment B  : " << HingeB0.position.transpose() << " - " << HingeB1.position.transpose()<<endl;
		cout<<"gab min " <<  gab_min<< "  Gab "<< Vec[k].Gab.transpose()<<endl;
		cout<<"r " <<  Vec[k].r.transpose() << "  r2 "<< Vec[k].r2.transpose()<<endl;
		switch (k)
			{
		case 0:
			cout<< " 0 segment-segment is shortest distance " <<endl;
			break;
		case 1:
			cout<< " 1 point-segment is shortest distance " <<endl;
			break;
		case 2:
			cout<< " 2 point-segment is shortest distance " <<endl;
			break;
		case 3:
			cout<< " 3 point-point is shortest distance " <<endl;
			break;
		case 4:
			cout<< " 4 point-point is shortest distance " <<endl;
			break;
		case 5:
			cout<< " 5 point-point is shortest distance " <<endl;
			break;
		case 6:
			cout<< " 6 point-point is shortest distance " <<endl;
		break;
			default:
				break;
			}
		cout<<endl;
	}
		

	


}






// TO DO Change type Vector3D to Type hinge // done
void computeExludedVolumeForces( Hinge& HingeA0,  const Hinge& HingeA1, Hinge& HingeB0, const Hinge& HingeB1, double threshold, double Excl_vol_fac, double bead_radius, double & gab_min, vector<Distance> & Vec, bool isForNeighborsList){  

	// TO DO, this is too verbose, make it more compact

	//vector<Distance> Vec(7);
	Distance tempVec;

	Distance_segments(HingeA0.position,HingeA1.position,HingeB0.position,HingeB1.position, tempVec.r, tempVec.r2, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	Vec[0]=tempVec; // 1

	Distance_point_segment(HingeA0.position, HingeB0.position, HingeB1.position, tempVec.r, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	tempVec.r2.setZero();
	Vec[1]=tempVec; // 2

	Distance_point_segment(HingeA1.position, HingeB0.position, HingeB1.position, tempVec.r, tempVec.Gab, tempVec.Gab_norm, tempVec.collisionCourse);
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	Vec[2]=tempVec; // 3

	Distance_point_point(HingeA0.position, HingeB0.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r.setZero();
	tempVec.r2.setZero();
	tempVec.collisionCourse =true; 
	Vec[3]=tempVec; // 4

	Distance_point_point(HingeA1.position, HingeB0.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r.setZero();
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	tempVec.collisionCourse =true;
	Vec[4]=tempVec; // 5

	Distance_point_point(HingeA0.position, HingeB1.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r = HingeB1.position-HingeB0.position;
	tempVec.r2.setZero();
	tempVec.collisionCourse =true;
	Vec[5]=tempVec; // 6

	Distance_point_point(HingeA1.position, HingeB1.position,tempVec.Gab, tempVec.Gab_norm);
	tempVec.r = HingeB1.position-HingeB0.position;
	tempVec.r2 =  HingeA1.position-HingeA0.position;
	tempVec.collisionCourse =true;
	Vec[6]=tempVec; // 7

	gab_min = Vec[6].Gab_norm;
	int k =6;

	for (int i = 0; i < 6; i++)
	{
		if (Vec[i].Gab_norm < gab_min  && Vec[i].collisionCourse)
		{
			gab_min = Vec[i].Gab_norm;
			k=i;

			


		}
	}
	/*
	if (isForNeighborsList)
	{
	switch (k)
			{
	case 0:
		cout<< " 0 segment-segment is shortest distance " <<endl;
		break;
	case 1:
		cout<< " 1 point-segment is shortest distance " <<endl;
		break;
	case 2:
		cout<< " 2 point-segment is shortest distance " <<endl;
		break;
	case 3:
		cout<< " 3 point-point is shortest distance " <<endl;
		break;
	case 4:
		cout<< " 4 point-point is shortest distance " <<endl;
		break;
	case 5:
		cout<< " 5 point-point is shortest distance " <<endl;
		break;
	case 6:
		cout<< " 6 point-point is shortest distance " <<endl;
		break;
			default:
				break;
			}
	}
	*/
	if (!isForNeighborsList && gab_min < threshold && Vec[k].collisionCourse )
	{
		Vector3d Excluded_volume_partial_force = -Vec[k].Gab * Excl_vol_fac * exp(-2 *(Vec[k].Gab_norm / bead_radius -2));
		HingeB0.exluded_volume_force +=  Excluded_volume_partial_force ;
		HingeA0.exluded_volume_force -=  Excluded_volume_partial_force ;
		HingeB0.excluded_volume_torque += Vec[k].r.cross(Excluded_volume_partial_force);
		HingeA0.excluded_volume_torque += Vec[k].r2.cross(-Excluded_volume_partial_force);
		/*
		cout<<" segment A  : " << HingeA0.position.transpose() << " - " << HingeA1.position.transpose()<<endl;
		cout<<" segment B  : " << HingeB0.position.transpose() << " - " << HingeB1.position.transpose()<<endl;
		cout<<" segment A Excluded volume force : "<< HingeA0.exluded_volume_force.transpose() << endl;
		cout<<" segment B Excluded volume force : "<< HingeB0.exluded_volume_force.transpose() << endl;
		*/
	}


}

double dotProduct(double Ax, double Ay, double Az, double Bx, double By, double Bz  ){

	return Ax*Bx + Ay*By + Az*Bz;

}

void findClosestPointSegmentDistance( double pointX, double pointY, double pointZ, double segmentAX, double segmentAY, double segmentAZ, double segmentBX, double segmentBY, double segmentBZ, double &distance, double &closestvectorX, double & closestvectorY){

	double ABX = segmentBX - segmentAX;
	double ABY = segmentBY - segmentAY;
	double ABZ = segmentBZ - segmentAZ;

	double t = dotProduct(  (pointX-segmentAX),(pointY-segmentAY),(pointZ-segmentAZ) , ABX , ABY, ABZ ) / dotProduct( ABX, ABY, ABZ , ABX, ABY, ABZ   );

	// clip it
	if (t<0){ t = 0;} 
	if (t>0){ t = 1;}

	double closestX = segmentAX + t * ABX;	  
	double closestY = segmentAY + t * ABY;
	double closestZ = segmentAZ + t * ABZ;

	closestvectorX = pointX-closestX;
	closestvectorY = pointY-closestY;


	 distance = sqrt( closestvectorX*closestvectorX + closestvectorY*closestvectorY    );
}


void calculateWallInteraction(Hinge & hinge1, double x1, double y1, double x2, double y2, double fiberRadius,  double excludedVolumeForce, double threshold){

	double closestDistance;
	double closestXVector;
	double closestYVector;


	findClosestPointSegmentDistance(hinge1.position[0],hinge1.position[1],0,x1,y1,0,x2,y2,0,closestDistance,closestXVector,closestYVector);
	//vector has correct sign, away from the wall
	if (closestDistance <= threshold)
	{
		

		hinge1.exluded_volume_force(0)+= closestXVector / closestDistance * excludedVolumeForce *  exp(-2 *(closestDistance / fiberRadius -2));
		hinge1.exluded_volume_force(1)+= closestYVector / closestDistance * excludedVolumeForce *  exp(-2 *(closestDistance / fiberRadius -2));
		//cout<<"Collision with wall, force in y is "<<hinge1.exluded_volume_force(1)<< "  closestDistance  "<<  closestDistance << " closestYVector "<<closestYVector << endl;
		//cout<< " x1  " <<x1<<" y1 " << y1<< " x2 "<<x2<<" y2 "<<y2<<" hinge x " <<hinge1.position[0]<< " hinge yx  " <<hinge1.position[1]<<endl;
	}


}



//collisions with walls
void interactionsWallsTri(vector<fiber> & fibers , int numNodes,vector<double> & meshPoints , vector<int> &meshConnectivity , vector<int> meshNeighbors , double fiberRadius , double excludedVolumeForce, double threshold, int nodesPerElement){
	#pragma omp parallel for
	for (int fiberID = 0; fiberID < fibers.size(); fiberID++)
	{
		for (int hingeID = 0; hingeID < fibers[fiberID].numberOfSegments; hingeID++)
		{
			for (int faceID = 0; faceID < 3; faceID++)// each face of the triangle
			{
				if (fibers[fiberID].hinges[hingeID].elementID!= -1 ) // the hinge is associated to an element
				{
					if (meshNeighbors[fibers[fiberID].hinges[hingeID].elementID*3+faceID ] ==-1){//-1 means that theres no neighbors in that face, thus is a boundary

						// the nodes that make that edge depndes on the index of the face:
						//face 0: nodes 0, 1, 
						//face 1: nodes 1, 2
						//face 2: nodes 2, 0
						int nodeA ;
							int nodeB;
						switch (faceID)
						{
						case 0:
							nodeA = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 0];
							nodeB = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 1];
							
							break;
						case 1:
							nodeA = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 1];
							nodeB = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 2];
								break;
							case 2:
								nodeA = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 2];
							nodeB = meshConnectivity[ fibers[fiberID].hinges[hingeID].elementID * nodesPerElement + 0];
								break;



						default:
							break;
						}

						calculateWallInteraction(fibers[fiberID].hinges[hingeID], meshPoints[ nodeA  ], meshPoints[ nodeA + numNodes  ],  meshPoints[ nodeB  ],meshPoints[ nodeB + numNodes  ],fiberRadius,excludedVolumeForce,threshold);
						

					}


				}
				



			}



		}
	}





}



void refineNeighbors(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  potentialCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices) {

	
	double threshold = 4  * beadRadius;
	double Excl_volume_face =1000000;
	vector<Distance> Vec(7); // used in compute exluded forces
	//cout<< " the threshold is " << threshold <<endl;
	// for every fiber
	int numberOFCollisionsMy= 0;
	int timeA;
	
	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
				 timeA =getMilliCount();
	#endif
	/*
	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		// for every segment of fiber i
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			
			// lets uppdate the neighbor list, naive programing, will be improved afterwards 
			// I dont know exactly how it will end up looking like
			vector<int> newNeighborFiberIndex;
			vector<int> newNeighborSegmentIndex;
			int newNumberOfNeighbours=0;
			for (unsigned int k = 0; k < fibers[i].segments[j].neighboringFibers.size(); k++)
			{

				newNeighborFiberIndex.resize(fibers[i].segments[j].neighboringFibers.size());
				newNeighborSegmentIndex.resize(fibers[i].segments[j].neighboringFibers.size());

				double min_Gab;
				int neighbouringFiberIndex =fibers[i].segments[j].neighboringFibers[k];
				int neighbouringSegmentIndex = fibers[i].segments[j].neighboringSegment[k];
				//TODO add flag for neighbors list
				computeExludedVolumeForces(fibers[i].hinges[j], fibers[i].hinges[j+1], fibers[neighbouringFiberIndex].hinges[neighbouringSegmentIndex], 
					fibers[neighbouringFiberIndex].hinges[neighbouringSegmentIndex+1], threshold, Excl_volume_face ,beadRadius, min_Gab, Vec,true);
				if (neighbouringFiberIndex != i && min_Gab < threshold  ) // != for non unique interactions, > for unique interactions
				{
					newNeighborFiberIndex[newNumberOfNeighbours] = neighbouringFiberIndex;
					newNeighborSegmentIndex[newNumberOfNeighbours] = neighbouringSegmentIndex;
					//newNeighborFiberIndex.push_back(neighbouringFiberIndex);
					//newNeighborSegmentIndex.push_back(neighbouringSegmentIndex);
					newNumberOfNeighbours++;
					numberOFCollisionsMy ++;

				}else
				{
					//cout<< " neighbor dismissed " <<endl;
				}
			}
			newNeighborFiberIndex.resize(newNumberOfNeighbours);
			newNeighborSegmentIndex.resize(newNumberOfNeighbours);

			//cout<< "1.  number of neighbours "<< fibers[i].segments[j].neighboringFibers.size() << " new number " << newNeighborFiberIndex.size()<<endl;
			fibers[i].segments[j].neighboringFibers.clear();
			fibers[i].segments[j].neighboringFibers = newNeighborFiberIndex;
			fibers[i].segments[j].neighboringSegment.clear();
			fibers[i].segments[j].neighboringSegment = newNeighborSegmentIndex;
			//cout<< "1.  number of neighbours "<< fibers[i].segments[j].neighboringFibers.size() << " new number " << newNeighborFiberIndex.size()<<endl;
		}
	}
	

	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
	cout << "time required to find neighbors my method " << getMilliSpan(timeA)<<" ms"<<endl;
	#endif
	*/

	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
				 timeA =getMilliCount();
	#endif

	thrust::host_vector<int> stencil;

	for (unsigned int i = 0; i < potentialCollisions.size(); i++)
	{
		int particleA = (long)(potentialCollisions[i]>>32); // global index of particle A
		int particleB = (long) potentialCollisions[i]; // global index of particle B
		int particleAfiberIndex = fiberIndices[particleA]; // fiber to which particle A belongs
		int particleBfiberIndex = fiberIndices[particleB]; // fiber to which particle B belongs
		int particleASegmentIndex = particleA - segmentIndices[particleAfiberIndex];// local index of particle A
		int particleBSegmentIndex = particleB - segmentIndices[particleBfiberIndex];// local index of particle B

		double min_Gab;
		computeExludedVolumeForces(fibers[particleAfiberIndex].hinges[particleASegmentIndex], fibers[particleAfiberIndex].hinges[particleASegmentIndex+1], 
						fibers[particleBfiberIndex].hinges[particleBSegmentIndex], fibers[particleBfiberIndex].hinges[particleBSegmentIndex+1],
						threshold, Excl_volume_face ,beadRadius, min_Gab, Vec,true);

		
		//cout<<"min_gab " << min_Gab << "threshold " <<threshold<<endl;
		if (min_Gab < threshold  )
				{
					
					stencil.push_back(1);
					//cout<<"do not delete interaction " <<endl;
				}
		else {
			stencil.push_back(0);
			//cout<<" delete interaction " <<endl;
		}

		 

	}
	Broadphase broadphaseManager;

	
	
	if (potentialCollisions.size() != 0 )
	{
		int numberOfCollisions = broadphaseManager.removeSameFiberCollisionsH(potentialCollisions,stencil);
		int realCollissions = stencil.size() - numberOfCollisions;
		potentialCollisions.resize(realCollissions);
	}
	#ifdef PRINT_DEBUG_INFO_NEIGHBORS_TIMING
	cout << "time required to find neighbors broadphase method " << getMilliSpan(timeA)<<" ms"<<endl;
	#endif

	
}




void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius) {

	
	double threshold = 10 * beadRadius;
	double Excl_volume_face =100;
	vector<Distance> Vec(7); // used in compute exluded forces
	//cout<< " the threshold is " << threshold <<endl;
	// for every fiber

	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		for (unsigned int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			fibers[i].hinges[j].exluded_volume_force.fill(0);
			fibers[i].hinges[j].torque.fill(0);
			fibers[i].hinges[j].excluded_volume_torque.fill(0);
		}
	}


	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		// for every segment of fiber i
		for (int j = 0; j < fibers[i].numberOfSegments; j++)
		{
			
			for (unsigned int k = 0; k < fibers[i].segments[j].neighboringFibers.size(); k++)
			{


				double min_Gab;
				int neighbouringFiberIndex =fibers[i].segments[j].neighboringFibers[k];
				int neighbouringSegmentIndex = fibers[i].segments[j].neighboringSegment[k];
				if (neighbouringFiberIndex != i) // If the neighbor does not belong to the same fiber as the current one
				{
				//TODO add flag for neighbors list // done
				computeExludedVolumeForces(fibers[i].hinges[j], fibers[i].hinges[j+1], fibers[neighbouringFiberIndex].hinges[neighbouringSegmentIndex], 
					fibers[neighbouringFiberIndex].hinges[neighbouringSegmentIndex+1], threshold, Excl_volume_face ,beadRadius, min_Gab,Vec,false);
				}
				
	
			}
			
		}
	}
}


void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  potentialCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices  ){

	double threshold = 2.1 * beadRadius;
	double Excl_volume_face =100;
	vector<Distance> Vec(7); // used in compute exluded forces
	//cout<< " the threshold is " << threshold <<endl;
	// for every fiber

	#ifdef PRINT_DEBUG_INFO
	cout<<" make sure that everything is zeroed "<<endl;
	#endif

	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			fibers[i].hinges[j].exluded_volume_force.fill(0);
			fibers[i].hinges[j].torque.fill(0);
			fibers[i].hinges[j].excluded_volume_torque.fill(0);

			#ifdef PRINT_DEBUG_INFO
			cout<<"fiber " << i << " segment " << j << endl;
			cout<< " excluded volume force " << fibers[i].hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fibers[i].hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fibers[i].hinges[j].excluded_volume_torque<<endl;
			#endif

		}
	}


	for (unsigned int i = 0; i < potentialCollisions.size(); i++)
	{
		int particleA = (long)(potentialCollisions[i]>>32); // global index of particle A
		int particleB = (long) potentialCollisions[i]; // global index of particle B
		int particleAfiberIndex = fiberIndices[particleA]; // fiber to which particle A belongs
		int particleBfiberIndex = fiberIndices[particleB]; // fiber to which particle B belongs
		int particleASegmentIndex = particleA - segmentIndices[particleAfiberIndex];// local index of particle A
		int particleBSegmentIndex = particleB - segmentIndices[particleBfiberIndex];// local index of particle B
		bool isForNeighborsList = false;
		double min_Gab;
		computeExludedVolumeForces(fibers[particleAfiberIndex].hinges[particleASegmentIndex], fibers[particleAfiberIndex].hinges[particleASegmentIndex+1], 
						fibers[particleBfiberIndex].hinges[particleBSegmentIndex], fibers[particleBfiberIndex].hinges[particleBSegmentIndex+1],
						threshold, Excl_volume_face ,beadRadius, min_Gab, Vec,isForNeighborsList);

		#ifdef PRINT_DEBUG_INFO
		cout<< "fiber " <<particleAfiberIndex <<" segment "<<particleASegmentIndex<<" force "<<	fibers[particleAfiberIndex].hinges[particleASegmentIndex].exluded_volume_force <<endl;
		cout<< "fiber " <<particleBfiberIndex <<" segment "<<particleBSegmentIndex<<" force "<<	fibers[particleBfiberIndex].hinges[particleBSegmentIndex].exluded_volume_force <<endl;

		#endif



	}



	#ifdef PRINT_DEBUG_INFO
	cout<<" results form calculations"<<endl;
	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			cout<<"fiber " << i << " segment " << j << endl;
			cout<< " excluded volume force " << fibers[i].hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fibers[i].hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fibers[i].hinges[j].excluded_volume_torque<<endl;
		}
	}
	#endif

}


void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  potentialCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices , int timeStep ){

	double threshold = 2.1 * beadRadius;

	

	double Excl_volume_face =100;
	vector<Distance> Vec(7); // used in compute exluded forces
	//cout<< " the threshold is " << threshold <<endl;
	// for every fiber

	#ifdef PRINT_DEBUG_INFO
	cout<<" make sure that everything is zeroed "<<endl;
	#endif

	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			fibers[i].hinges[j].exluded_volume_force.fill(0);
			fibers[i].hinges[j].torque.fill(0);
			fibers[i].hinges[j].excluded_volume_torque.fill(0);

			#ifdef PRINT_DEBUG_INFO
			cout<<"fiber " << i << " segment " << j << endl;
			cout<< " excluded volume force " << fibers[i].hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fibers[i].hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fibers[i].hinges[j].excluded_volume_torque<<endl;
			#endif

		}
	}


	for (unsigned int i = 0; i < potentialCollisions.size(); i++)
	{
		int particleA = (long)(potentialCollisions[i]>>32); // global index of particle A
		int particleB = (long) potentialCollisions[i]; // global index of particle B
		int particleAfiberIndex = fiberIndices[particleA]; // fiber to which particle A belongs
		int particleBfiberIndex = fiberIndices[particleB]; // fiber to which particle B belongs
		int particleASegmentIndex = particleA - segmentIndices[particleAfiberIndex];// local index of particle A
		int particleBSegmentIndex = particleB - segmentIndices[particleBfiberIndex];// local index of particle B
		bool isForNeighborsList = false;
		double min_Gab;
		computeExludedVolumeForces(fibers[particleAfiberIndex].hinges[particleASegmentIndex], fibers[particleAfiberIndex].hinges[particleASegmentIndex+1], 
						fibers[particleBfiberIndex].hinges[particleBSegmentIndex], fibers[particleBfiberIndex].hinges[particleBSegmentIndex+1],
						threshold, Excl_volume_face ,beadRadius, min_Gab, Vec,isForNeighborsList,timeStep);

		#ifdef PRINT_DEBUG_INFO
		cout<< "fiber " <<particleAfiberIndex <<" segment "<<particleASegmentIndex<<" force "<<	fibers[particleAfiberIndex].hinges[particleASegmentIndex].exluded_volume_force <<endl;
		cout<< "fiber " <<particleBfiberIndex <<" segment "<<particleBSegmentIndex<<" force "<<	fibers[particleBfiberIndex].hinges[particleBSegmentIndex].exluded_volume_force <<endl;

		#endif



	}



	#ifdef PRINT_DEBUG_INFO
	cout<<" results form calculations"<<endl;
	for (unsigned int i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			cout<<"fiber " << i << " segment " << j << endl;
			cout<< " excluded volume force " << fibers[i].hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fibers[i].hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fibers[i].hinges[j].excluded_volume_torque<<endl;
		}
	}
	#endif

}