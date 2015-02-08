
#include "bending torque.h"



Vector3d bendingTorque( Vector3d point1, Vector3d point2, Vector3d point3, double E_young, double Inertia_moment){
	// computes the torque from two adjacent rods, or three adjacent hinges
	
	double angle;
	Vector3d T_b;
	Vector3d vector1 = point2 - point1;
	Vector3d vector2 = point3 - point2;

	double length1 = sqrt(abs(vector1.dot(vector1)));
	double length2 = sqrt(abs(vector2.dot(vector2)));

	

	double K_bending = (2 * E_young * Inertia_moment) / (length1+length2);
	
	double angleCos = vector1.dot(vector2)/(length1*length2);
	
	if (angleCos >= 1 )
	{
		angle = 0;
	} else if (angleCos <=-1 )
	{
		angle = pi/2;
	} else
	{
		angle =acos(angleCos);
	}

	
	Vector3d crossProduct = vector1.cross(vector2);
	
	Vector3d e_theta = crossProduct /(sqrt(crossProduct.dot(crossProduct)));
	
	if ( abs(crossProduct.x()) <= DBL_EPSILON && 
		abs(crossProduct.y())  <= DBL_EPSILON &&
		abs(crossProduct.z())  <= DBL_EPSILON){
			T_b.fill(0);
	}else{
		T_b = K_bending *angle*e_theta;
	}
	/*
	cout<< " vector1 : " << vector1<<endl;
	cout<< " vector2 : " << vector2<<endl;
	cout<< "k_bending " << K_bending << endl;
	cout<< "angleCos " << angleCos << endl;
	cout<< "angle " << angle << endl;
	cout<< "crossProduct " << crossProduct << endl;
	cout<< "e_theta " << e_theta << endl;
	cout<< "angleCos " << angleCos << endl;
	*/
	return T_b;

}



void bendingTorqueForAllFibers(vector<fiber> & fibers, vector<Hinge> & hinges,double E_young, double inertia_moment ){

	// copying literally Daniels routines




	for (int i = 0; i < fibers.size(); i++)
	{
		for (int j = 0; j < fibers[i].numberOfHinges; j++)
		{
			fibers[i].hinges[j].torque.fill(0);
			fibers[i].hinges[j].angle =0;
		}
	}


	for (int i = 0; i < fibers.size(); i++)
	{

		//cout<< " fiber ID " << i << endl;

		if (fibers[i].numberOfHinges >=3)
		{
			// Using old data structure
			/*
			for (int j = fibers[i].firstHinge+1; j < fibers[i].firstHinge+fibers[i].numberOfHinges-2 ; j++)
			{
				
				Vector3d T_b = bendingTorque(hinges[j-1].position , hinges[j].position , hinges[j+1].position, E_young , inertia_moment );
				hinges[j-1].torque += T_b;
				hinges[j].torque -=T_b;

				//cout<< " Torque old data structure " << T_b.transpose() << endl; 

			}
			*/
				//cout<<"old "<<fibers[i].firstHinge+fibers[i].numberOfHinges-2 -(fibers[i].firstHinge+1)<<endl;

			// Using new data Structure
			for (int j = 1; j < fibers[i].numberOfHinges-1; j++) // -1 because we are skipping the first and last hinge
			{
				Vector3d T_b = bendingTorque(fibers[i].hinges[j-1].position , 
											 fibers[i].hinges[j].position, 
											 fibers[i].hinges[j+1].position, 
											 E_young, inertia_moment);
				//cout<<" torque " << T_b <<endl;
				fibers[i].hinges[j-1].torque += T_b;
				fibers[i].hinges[j].torque -= T_b;
				//cout<< " Torque new data structure " << T_b.transpose() << endl; 
			}

			//cout<<"new "<<fibers[i].numberOfHinges-2-1 << endl;
			 

		}
	}
}