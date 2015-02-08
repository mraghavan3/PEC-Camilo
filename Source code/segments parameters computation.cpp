#include "SEgments parameters computation.h"


void compute_velocities_mesh(Vector3d position, Vector3d &velocity, Vector3d & vorticity, Hinge hinge0, Hinge hinge1 ){

	
	double distance1 = sqrt ( (hinge0.position[0] - hinge1.position[0]) * (hinge0.position[0] - hinge1.position[0]) + (hinge0.position[1] - hinge1.position[1]) *(hinge0.position[1] - hinge1.position[1]) + (hinge0.position[2] - hinge1.position[2]) *(hinge0.position[2] - hinge1.position[2]));
	double distance2 = sqrt ( (position[0] - hinge1.position[0]) * (position[0] - hinge1.position[0]) + (position[1] - hinge1.position[1]) *(position[1] - hinge1.position[1]) + (position[2] - hinge1.position[2]) *(position[2] - hinge1.position[2]));

	double weight = distance2-distance1; 
	double weight2 = 1-weight;
	velocity(0) = weight*hinge1.velocity[0] + weight2*hinge0.velocity[0];
	velocity(1) = weight*hinge1.velocity[1] + weight2*hinge0.velocity[1];
	velocity(2) = weight*hinge1.velocity[2] + weight2*hinge0.velocity[2];
	vorticity(0) = weight*hinge1.vorticity[0] + weight2*hinge0.vorticity[0];
	vorticity(1) = weight*hinge1.vorticity[1] + weight2*hinge0.vorticity[1];
	vorticity(2) = weight*hinge1.vorticity[2] + weight2*hinge0.vorticity[2];

	//cout<<" compute velocity for a mesh x " << velocity(0) << " y " << velocity(1) <<endl;
}

void compute_velocities_mesh2(double distance1, double distance2, Vector3d &velocity, Vector3d & vorticity, Hinge hinge0, Hinge hinge1 ){

	
	double weight = distance2/distance1; 
	double weight2 = 1-weight;
	velocity(0) = weight*hinge1.velocity[0] + weight2*hinge0.velocity[0];
	velocity(1) = weight*hinge1.velocity[1] + weight2*hinge0.velocity[1];
	velocity(2) = weight*hinge1.velocity[2] + weight2*hinge0.velocity[2];
	vorticity(0) = weight*hinge1.vorticity[0] + weight2*hinge0.vorticity[0];
	vorticity(1) = weight*hinge1.vorticity[1] + weight2*hinge0.vorticity[1];
	vorticity(2) = weight*hinge1.vorticity[2] + weight2*hinge0.vorticity[2];

	//cout<<" compute velocity for a mesh x " << velocity(0) << " y " << velocity(1) <<endl;
}

void compute_velocities(Vector3d position, Vector3d &velocity, Vector3d &omega ){

	// To Do, appropiate interpolation
	// simple shear flow analytical solution
	// velocity components in y and z are 0
	velocity(1) = 0;
	velocity(2) = 0;

	// shear rate = macVel / ( Ymax - Ymin)
	double maxVelocity = 0.2;
	double Ymax = 0.001999;
	double Ymin = 0;
	double shearRate = 0.5;

	velocity(0) = shearRate * position.y();
	//cout<<"velocity "<<velocity(1) << " position Y" << position.y();
	omega(0) =0;
	omega(1) = 0;
	omega(2) = -0.5*shearRate; //2D flow vorticity is in the Z direction

}


void hinge_parameter_computation(Hinge & hinge1, Hinge hinge2, double beadRadius){


	hinge1.fluid_velocity_sum.fill(0);
	hinge1.sumOmegaFluid.fill(0);
	hinge1.r_sum.fill(0);
	hinge1.r_product_sum.fill(0);
	hinge1.r_times_u_sum.fill(0);

	
	hinge1.r = hinge2.position - hinge1.position;
	double segmentLength = sqrt(hinge1.r.dot(hinge1.r));
	Vector3d unitaryVector = hinge1.r / segmentLength;
	hinge1.numberOfBeads = 2*floor((segmentLength/2)/(2*beadRadius));
	Vector3d segmentMiddle = (hinge1.position + hinge2.position)/2;

	Vector3d startingPosition = segmentMiddle - unitaryVector * (2*beadRadius * hinge1.numberOfBeads/2-beadRadius);
	//cout<<startingPosition<<endl;
	// I dont know what this means, I took it directly from Daniels code
	Vector3d X_j;
	Vector3d X_local;
	Vector3d X_local_matrix;
	Vector3d X_Local_matrix_T; // T is for transposed
	Vector3d velocity;
	Vector3d omega;
	//cout<< " number of beads " << hinge1.numberOfBeads << " hinge r " << hinge1.r << " segment length " << segmentLength;


	for (int i = 0; i < hinge1.numberOfBeads; i++)
	{
		
		X_j = startingPosition + (i)*unitaryVector*(2*beadRadius);
		X_local = X_j - hinge1.position;

		
		compute_velocities (X_j , velocity , omega );
		hinge1.fluid_velocity_sum += velocity;
		hinge1.sumOmegaFluid += omega;
		hinge1.r_sum += X_local;

		hinge1.r_product_sum +=  X_local * X_local.transpose() ;
		hinge1.r_times_u_sum += X_local *  velocity.transpose();//( velocity* X_local.transpose()) ;//
	}
	/*
	cout<< " X_local : " << X_local << endl;
	cout<< " velocity.transpose() : " << velocity.transpose() << endl;
	*/
	hinge1.averageviscosity = 11.6; // To do change this to the interpolation scheme
	hinge1.torque.fill(0);
}



void hinge_parameter_computation_mesh(Hinge & hinge1, Hinge hinge2, double beadRadius){


	hinge1.fluid_velocity_sum.fill(0);
	hinge1.sumOmegaFluid.fill(0);
	hinge1.r_sum.fill(0);
	hinge1.r_product_sum.fill(0);
	hinge1.r_times_u_sum.fill(0);

	
	hinge1.r = hinge2.position - hinge1.position;
	double segmentLength = sqrt(hinge1.r.dot(hinge1.r));
	Vector3d unitaryVector = hinge1.r / segmentLength;
	hinge1.numberOfBeads = 2*floor((segmentLength/2)/(2*beadRadius));
	Vector3d segmentMiddle = (hinge1.position + hinge2.position)/2;

	Vector3d startingPosition = segmentMiddle - unitaryVector * (2*beadRadius * hinge1.numberOfBeads/2-beadRadius);
	//cout<<startingPosition<<endl;
	// I dont know what this means, I took it directly from Daniels code
	Vector3d X_j;
	Vector3d X_local;
	Vector3d X_local_matrix;
	Vector3d X_Local_matrix_T; // T is for transposed
	Vector3d velocity;
	Vector3d omega;
	//cout<< " number of beads " << hinge1.numberOfBeads << " hinge r " << hinge1.r << " segment length " << segmentLength;


	for (int i = 0; i < hinge1.numberOfBeads; i++)
	{
		
		X_j = startingPosition + (i)*unitaryVector*(2*beadRadius);
		X_local = X_j - hinge1.position;

		
		//compute_velocities_mesh(X_j,velocity,omega,hinge1,hinge2);
		
		compute_velocities_mesh2(segmentLength,i*2*beadRadius,velocity,omega,hinge1,hinge2);
		hinge1.fluid_velocity_sum += velocity;
		hinge1.sumOmegaFluid += omega;
		hinge1.r_sum += X_local;

		hinge1.r_product_sum +=  X_local * X_local.transpose() ;
		hinge1.r_times_u_sum += X_local *  velocity.transpose();//( velocity* X_local.transpose()) ;//
	}
	/*
	cout<< " X_local : " << X_local << endl;
	cout<< " velocity.transpose() : " << velocity.transpose() << endl;
	*/
	hinge1.averageviscosity = 100; // To do change this to the interpolation scheme
	hinge1.torque.fill(0);
}


void fibers_parameters_computation(vector<fiber> & fibers, double beadRadius){

	#pragma omp parallel for
	for (int i = 0; i < fibers.size(); i++) // for each fiber
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++) // for each hinge in the fiber
		{
			hinge_parameter_computation(fibers[i].hinges[j], fibers[i].hinges[j+1], beadRadius);
			
		}

	}



}


void interpolate3Tri(double x, double y, double z, double & fx, double & fy, double & fz,
					double x1, double y1, double z1, double fx1, double fy1, double fz1,
					double x2, double y2, double z2, double fx2, double fy2, double fz2,
					double x3, double y3, double z3, double fx3, double fy3, double fz3
					){

// Function returns the interpolated vector function (fx,fy,fz) evaluated in a point (x,y,z), 
// using a linear Triangle
// the function can be interpolated with:
//			ie. fx = Xi1*fx1 +  Xi2*fx2 + Xi3*fx3 

double x12 = x1-x2; double x13 = x1-x3; 
double x21 = x2-x1 ; double x23 = x2-x3; 
double x31 = x3-x1; double x32 = x3-x2; 

double y12 = y1-y2; double y13 = y1-y3; 
double y21 = y2-y1 ; double y23 = y2-y3; 
double y31 = y3-y1; double y32 = y3-y2; 
/*
double z12 = z1-z2; double z13 = z1-z3; 
double z21 = z2-z1 ; double z23 = z2-z3; 
double z31 = z3-z1; double z32 = z3-z2; 
*/
double A = ((x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1));

//cout<< " A " <<A <<endl;


double A23 = x2*y3-x3*y2;
double A31 = x3*y1-x1*y3;
double A12 = x1*y2-x2*y1;

double Xi1 = (A23+ y23*x + x32*y  ) /A  ; 
double Xi2 = (A31 +y31*x + x13*y  ) /A  ;
double Xi3 = (A12 + y12*x + x21*y  ) /A  ; 

//cout<< " x1 " <<x1<< " x2 " <<x2<< " x3 " <<x3 <<endl;
//cout<< " y1 " <<y1<< " y2 " <<y2<< " y3 " <<y3 <<endl;
//cout<< " z1 " <<z1<< " z2 " <<z2<< " z3 " <<z3 <<endl;

//cout<< " Xi1 " <<Xi1<< " Xi2 " <<Xi2<< " Xi3 " <<Xi3 <<endl;

fx = Xi1*fx1 +  Xi2*fx2 + Xi3*fx3 ;
fy = Xi1*fy1 +  Xi2*fy2 + Xi3*fy3 ;
fz = Xi1*fz1 +  Xi2*fz2 + Xi3*fz3 ;

/*
double a1 =x2*y3-y2*x3;
double a2 =x3*y1-y3*x1;
double a3 =x1*y2-y1*x2;

double A = a1 + a2 +a3;

a1 = a1/(2*A);
a2 = a2/(2*A);
a3 = a3/(2*A);

double b1  = (y2-y3)/(2*A);
double b2  = (y3-y1)/(2*A);
double b3  = (y1-y2)/(2*A);

double c1  = (x3-x2)/(2*A);
double c2  = (x1-x3)/(2*A);
double c3  = (x2-x1)/(2*A);

double n1 = a1 + b1*x + c1*y;
double n2 = a2 + b2*x + c2*y;
double n3 = a3 + b3*x + c3*y;

fx = n1*fx1 +  n2*fx2 + n3*fx3 ;
fy = n1*fy1 +  n2*fy2 + n3*fy3 ;
fz = n1*fz1 +  n2*fz2 + n3*fz3 ;

*/
}

void interpolate6Tri(double &x, double &y, double &z, double & fx, double & fy, double & fz,
					double &x1, double &y1, double& z1, double& fx1, double &fy1, double &fz1,
					double& x2, double &y2, double &z2, double &fx2, double &fy2, double &fz2,
					double& x3, double& y3, double &z3, double &fx3, double& fy3, double &fz3,
					double &x4, double& y4, double& z4, double &fx4, double &fy4, double &fz4,
					double& x5, double &y5, double &z5, double &fx5, double &fy5, double& fz5,
					double &x6, double &y6, double& z6, double &fx6, double &fy6, double& fz6
					){

// Function returns the interpolated vector function (fx,fy,fz) evaluated in a point (x,y,z), 
// using a linear Triangle
// the function can be interpolated with:
//			ie. fx = N1*fx1 +  N2*fx2 + N3*fx3 + N4*fx4, + N5*fx5 + N6*fx6 ;

/*double x12 = x1-x2;*/ double x13 = x1-x3; 
double x21 = x2-x1 ; //double x23 = x2-x3; 
/*double x31 = x3-x1;*/ double x32 = x3-x2; 

double y12 = y1-y2; //double y13 = y1-y3; 
/*double y21 = y2-y1 ;*/ double y23 = y2-y3; 
double y31 = y3-y1; //double y32 = y3-y2; 

// a tri is in 2D

/*
double z12 = z1-z2; double z13 = z1-z3; 
double z21 = z2-z1 ; double z23 = z2-z3; 
double z31 = z3-z1; double z32 = z3-z2; 
*/
double A = ((x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1))/2;

double A23 = x2*y3-x3*y2;
double A31 = x3*y1-x1*y3;
double A12 = x1*y2-x2*y1;

double Xi1 = (A23 + y23*x + x32*y  ) /A  ; 
double Xi2 = (A31 + y31*x + x13*y  ) /A  ;
double Xi3 = (A12 + y12*x + x21*y  ) /A  ; 

// a small numerical error should not be a big issue here
//if(abs(Xi1)<0.0001) {Xi2-=abs(Xi1);  Xi1=0; } 
//if(abs(Xi2)<0.0001) {Xi1-=abs(Xi2);Xi2=0; } 
//if(abs(Xi3)<0.0001) {Xi1-=abs(Xi3);Xi3=0; } 

double N1 = Xi1*(2*Xi1-1);
double N2 = Xi2*(2*Xi2-1);
double N3 = Xi3*(2*Xi3-1);

double N4 = 4*Xi1*Xi2;
double N5 = 4*Xi2*Xi3;
double N6 = 4*Xi3*Xi1;



fx = N1*fx1 +  N2*fx2 + N3*fx3 + N4*fx4, + N5*fx5 + N6*fx6 ;
fy = N1*fy1 +  N2*fy2 + N3*fy3 + N4*fy4, + N5*fy5 + N6*fy6 ;
fz = N1*fz1 +  N2*fz2 + N3*fz3 + N4*fz4, + N5*fz5 + N6*fz6 ;

}



void interpolateFunction(Hinge &hinge0, vector<double > &meshPoints,  vector<int> &meshConnectivity,  vector<double> &meshVelocities, vector<double>  meshVorticity,int numElements, int numNodes, int nodesPerElement  ){

	switch (nodesPerElement)
	{

	case 3: // linear elements
		//interpolate velocities
		/*cout<< " Element  ID " <<  hinge0.elementID<<endl;
		cout<< " Node 1 " << meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] ;
		cout<< " Node 2 " << meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] ;
		cout<< " Node 3 " << meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] <<endl;*/
		if (!hinge0.isStationary)
		{

		
		interpolate3Tri(hinge0.position[0], hinge0.position[1], hinge0.position[2], hinge0.velocity[0], hinge0.velocity[1],hinge0.velocity[2],
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// coords for node 0
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// vels for node 0
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// coords for node 1
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// vels for node 1
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// coords for node 2
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ]);// vels for node 2
		/*
		cout<< " element ID " << hinge0.elementID<<endl; 
		cout<<"velocity in x first node "<<meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ]<< " is " << meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ]<<endl;
		cout<<" velocity in x second node "<<meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ]<< " is "<< meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ]<<endl;
		cout<<" velocity in x third node "<<meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ]<< " is "<< meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ]<<endl;
		cout<<endl;
		cout<<"hinge velocity in x " << hinge0.velocity[0] << " hinge velocity in y " << hinge0.velocity[1] <<endl;*/
		//interpolate vorticities
		interpolate3Tri(hinge0.position[0], hinge0.position[1], hinge0.position[2], hinge0.vorticity[0], hinge0.vorticity[1],hinge0.vorticity[2],
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// coords for node 0
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// vels for node 0
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// coords for node 1
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// vels for node 1
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// coords for node 2
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ]);// vels for node 2
		
		}else{
			hinge0.velocity[0] = 0;
			hinge0.velocity[1] = 0;
			hinge0.velocity[2] = 0;

			hinge0.vorticity[0] = 0;
			hinge0.vorticity[1] = 0;
			hinge0.vorticity[2] = 0;


		
		}


		break;
	case 6: // quadratic elements
		//interpolate velocities
	interpolate6Tri(hinge0.position[0], hinge0.position[1], hinge0.position[2], hinge0.velocity[0], hinge0.velocity[1],hinge0.velocity[2],
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// coords for node 0
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// vels for node 0
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// coords for node 1
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// vels for node 1
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// coords for node 2
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// vels for node 2
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*2 ],// coords for node 3
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*2 ],// vels for node 3
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*2 ],// coords for node 4
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*2 ],// vels for node 4
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*2 ],// coords for node 5
		meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*0 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*1 ], meshVelocities[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*2 ]);// vels for node 5
	
	//interpolate vorticities

	interpolate6Tri(hinge0.position[0], hinge0.position[1], hinge0.position[2], hinge0.vorticity[0], hinge0.vorticity[1],hinge0.vorticity[2],
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// coords for node 0
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  0 ] + numNodes*2 ],// vels for node 0
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// coords for node 1
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  1 ] + numNodes*2 ],// vels for node 1
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// coords for node 2
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  2 ] + numNodes*2 ],// vels for node 2
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*2 ],// coords for node 3
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  3 ] + numNodes*2 ],// vels for node 3
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*2 ],// coords for node 4
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  4 ] + numNodes*2 ],// vels for node 4
		meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*0 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*1 ], meshPoints[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*2 ],// coords for node 5
		meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*0 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*1 ], meshVorticity[  meshConnectivity[ hinge0.elementID * nodesPerElement +  5 ] + numNodes*2 ]);// vels for node 5
	
	//cout<<"hinge velocity in x " << hinge0.velocity[0] << " hinge velocity in y " << hinge0.velocity[1] <<endl;


		break;
	default:
		break;
	}



	

}


void fibers_parameters_computation_mesh(vector<fiber> & fibers, double beadRadius, vector<double > &meshPoints, vector<int> &meshConnectivity, vector<double> &meshVelocities, vector<double> &meshVorticity, int numElements, int numNodes, int nodesPerElement ){
 
//#pragma omp parallel for
	for (int i = 0; i < fibers.size(); i++) // for each fiber
	{
		for (int j = 0; j < fibers[i].numberOfSegments; j++) // for each hinge in the fiber
		{
			if (j==0)
			{

				//interpolates velocities and vorticities using the mesh information, assigns value to each hinge
				interpolateFunction(fibers[i].hinges[j],meshPoints,meshConnectivity,meshVelocities,meshVorticity,numElements,numNodes,nodesPerElement);
				interpolateFunction(fibers[i].hinges[j+1],meshPoints,meshConnectivity,meshVelocities,meshVorticity,numElements,numNodes,nodesPerElement);
				/*cout<<"fiber " << i  << " hinge " << j ;
				cout<<" velocity " <<  fibers[i].hinges[j].velocity.transpose()<<endl;
				*///cout<<  " vorticity " <<  fibers[i].hinges[j].vorticity<<endl;
			}else
			{
				interpolateFunction(fibers[i].hinges[j+1],meshPoints,meshConnectivity,meshVelocities,meshVorticity,numElements,numNodes,nodesPerElement);
				/*cout<<"fiber " << i  << " hinge " << j ;
				cout<<" velocity " <<  fibers[i].hinges[j].velocity.transpose()<<endl;*/
			}

			hinge_parameter_computation_mesh(fibers[i].hinges[j], fibers[i].hinges[j+1], beadRadius);
			
		}

	}


}