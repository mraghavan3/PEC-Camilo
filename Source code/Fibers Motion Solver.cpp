#include "Fibers Motion Solver.h"
#include <omp.h>
// Use this function carefully, TO DO check bounds and assert things. 

void copyMatrixRange(MatrixXd & targetMatrix, int initRowTarget, int finalRowTarget,
					 int initColumnTarget, int finalColumnTarget,
					 MatrixXd sourceMatrix, int initRowSource, int initColumnSource){

						 //cout<< "number of rows to be copied " << finalRowTarget-initRowTarget+1<<endl;
						 //cout<< " number of columns to be copied " << finalColumnTarget - initColumnTarget+1<<endl;

						 for (int i = 0; i < finalRowTarget-initRowTarget+1; i++)
						 {
							 for (int j = 0; j < finalColumnTarget - initColumnTarget +1 ; j++)
							 {
								 targetMatrix(initRowTarget+i,initColumnTarget+j)=sourceMatrix(initRowSource+i,initColumnSource+j);

							 }
						 }
}



// copies vector range to another target
void copyVectorRange(VectorXd & targetVector, int initialIndex, int finalIndex, VectorXd sourceVector, int initialIndexSource){

	for (int i = 0; i < finalIndex-initialIndex+1; i++)
	{
		targetVector(initialIndex+i) = sourceVector(initialIndexSource+i);

	}
}

// for debugging
void miniMatrixAssembly(MatrixXd & miniMatrix, VectorXd &miniVector, const Hinge & hinge1, Hinge hinge2, double dragCoefficientVelocity, double dragCoefficientOmega, int timeStep){


	// when comparing to fortran code remember that indexing is 0 based here as oposed to fortran
	miniMatrix.fill(0);
	miniVector.fill(0);

	////////////////////
	miniMatrix(0,3) = 1;
	miniMatrix(0,7) = hinge1.r(2); // this is in Z
	miniMatrix(0,8) = -hinge1.r(1); // this is in Y
	miniMatrix(0,12) = -1;
	// I checked the numbers
	

	////////////////////
	miniMatrix(1,4) = 1;
	miniMatrix(1,6) = -hinge1.r(2); // this is in Z
	miniMatrix(1,8) = hinge1.r(0); // this is in X
	miniMatrix(1,13) = -1;
	// I checked the numbers

	
	////////////////////
	miniMatrix(2,5) = 1;
	miniMatrix(2,6) = hinge1.r(1); // this is in Z
	miniMatrix(2,7) = -hinge1.r(0); // this is in X
	miniMatrix(2,14) = -1;
	// I checked the numbers

	

	///////////////////////
	miniMatrix(3, 0)=  1;
	miniMatrix(3, 3)= -dragCoefficientVelocity* hinge1.averageviscosity * hinge1.numberOfBeads  ;
	miniMatrix(3, 7)= -dragCoefficientVelocity* hinge1.averageviscosity* hinge1.r_sum(2);
	miniMatrix(3, 8)= +dragCoefficientVelocity *hinge1.averageviscosity* hinge1.r_sum(1);
	miniMatrix(3,9)= -1;

	// I checked the numbers


	/******************/
	miniVector(3)= -dragCoefficientVelocity * hinge1.averageviscosity *  hinge1.fluid_velocity_sum(0) - hinge1.exluded_volume_force(0); // add excluded volume force
	// the number it returns is -0, I dont know if there might be any numerical implication to this. 

	

	/*
	cout<< " minimatrix (4,1) " << miniMatrix(3,0) << " (4,4) " <<miniMatrix(3,3);
	cout<< " (4,8) " << miniMatrix(3,7 ) << " (4,9)" << miniMatrix(3,8);
	cout<< " (4,10) " << miniMatrix(3,9 )<<endl;
	*/
	///////////////////

	miniMatrix(4, 1)=  1;
	miniMatrix(4, 4)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.numberOfBeads;
	miniMatrix(4, 6)= +dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(2);
	miniMatrix(4, 8)= -dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(0);
	miniMatrix(4,10)= -1;
	// correct result
	/*
	cout<< " minimatrix (5,2) " << miniMatrix(4,1) << " (5,5) " <<miniMatrix(4,4);
	cout<< " (5,7) " << miniMatrix(4,6 ) << " (5,9)" << miniMatrix(4,8);
	cout<< " (5,11) " << miniMatrix(4,10 )<<endl;
	*/
	//////////////////

	miniVector(4)= -dragCoefficientVelocity * hinge1.averageviscosity *  hinge1.fluid_velocity_sum(1)- hinge1.exluded_volume_force(1); // add excluded volume force
	// result is correct
	///////////////////

	miniMatrix(5, 2)=  1;
	miniMatrix(5, 5)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.numberOfBeads;
	miniMatrix(5, 6)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.r_sum(1);
	miniMatrix(5, 7)= +dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(0);
	miniMatrix(5,11)= -1;
	// result is correct
	/*
	cout<< " minimatrix (6,3) " << miniMatrix(5,2) << " (6,6) " <<miniMatrix(5,5);
	cout<< " (6,7) " << miniMatrix(5,6 ) << " (6,8)" << miniMatrix(5,7);
	cout<< " (6,12) " << miniMatrix(5,11 )<<endl;
	*/
	///////////////////////

	miniVector(5)= -dragCoefficientVelocity* hinge1.averageviscosity *  hinge1.fluid_velocity_sum(2)- hinge1.exluded_volume_force(2); // add excluded volume force

#ifdef PRINT_DEBUG_INFO
	cout<< "minivector(5) "<<miniVector(5)<<endl;
	cout<< "dragCoefficientVelocity "<<dragCoefficientVelocity<<endl;
	cout<< "hinge1.averageviscosity "<<hinge1.averageviscosity<<endl;
	cout<< "hinge1.fluid_velocity_sum(2) "<<hinge1.fluid_velocity_sum(2)<<endl;
	cout<< "hinge1.excludedVolumeForce(2) "<<hinge1.exluded_volume_force(2)<<endl;
	cout<< "hinge1.excludedVolumeForce "<<hinge1.exluded_volume_force<<endl;
#endif


	// again -0 as the answer
	
	///////////////////////

	miniMatrix(6,6)=1; // again!!!??? (6,6) is repeated
	miniMatrix(6, 4)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(2);
	miniMatrix(6, 5)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(1);
	miniMatrix(6, 6)= -dragCoefficientOmega * hinge1.averageviscosity *  hinge1.numberOfBeads  
					- dragCoefficientVelocity *hinge1.averageviscosity  *  ( hinge1.r_product_sum(2,2)+hinge1.r_product_sum(1,1)) ;           
	miniMatrix(6, 7)=  dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,1);
	miniMatrix(6, 8)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,2);
	miniMatrix(6,10)= +hinge1.r(2);
	miniMatrix(6,11)= -hinge1.r(1);
	// results are correct
	/*
	cout<< " minimatrix (7,7) " << miniMatrix(6,6) << " (7,5) " <<miniMatrix(6,4);
	cout<< " (7,6) " << miniMatrix(6,5 ) << " (7,7)" << miniMatrix(6,6);
	cout<< " (7,8) " << miniMatrix(6,7 ) << " (7,9)" << miniMatrix(6,8);
	cout<< " (7,11) " << miniMatrix(6,10 ) << " (7,12)" << miniMatrix(6,11)<<endl;
	*/
	/////////////////////////
	
	miniVector(6)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(0) 
		+dragCoefficientVelocity * hinge1.averageviscosity  * ( hinge1.r_times_u_sum(2,1)- hinge1.r_times_u_sum(1,2))
		-hinge1.torque(0) - hinge1.excluded_volume_torque(0);  // add excluded volume force 
	// correct result

	//////////////////////////

	miniMatrix(7, 3)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(2);
	miniMatrix(7, 5)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(0);
	miniMatrix(7, 6)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,1);
	miniMatrix(7, 7)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.numberOfBeads
           -dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_product_sum(2,2)+hinge1.r_product_sum(0,0));
	miniMatrix(7, 8)=  dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(1,2);
	miniMatrix(7,9)= -hinge1.r(2);
	miniMatrix(7,11)= +hinge1.r(0);
	// correct result 
	/*
	cout<< " minimatrix (8,4) " << miniMatrix(7,3) << " (8,6) " <<miniMatrix(7,5);
	cout<< " (8,7) " << miniMatrix(7,6 ) << " (8,8)" << miniMatrix(7,7);
	cout<< " (8,9) " << miniMatrix(7,8 ) << " (8,10)" << miniMatrix(7,9);
	cout<< " (8,12) " << miniMatrix(7,11)<<endl;
	*/
	///////////////////////////

	miniVector(7)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(1)
		+dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_times_u_sum(0,2)-hinge1.r_times_u_sum(2,0))
           -hinge1.torque(1)- hinge1.excluded_volume_torque(0);// add excluded volume torque 
	/*
	cout<<" dragCoefficientOmega " << dragCoefficientOmega;
	cout<<" averageviscosity " << hinge1.averageviscosity;
	cout<<" sumOmegaFluid(1) " << hinge1.sumOmegaFluid(1);
	cout<<" r_times_u_sum(0,2) " << hinge1.r_times_u_sum(0,2);
	cout<<" r_times_u_sum(2,0) " << hinge1.r_times_u_sum(2,0);
	cout<<" minivector(7) " << miniVector(7);
	cout<<endl;
	*/
	// correct result
	
	//////////////////////////////

	miniMatrix(8, 3)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(1);
	miniMatrix(8, 4)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(0);
	miniMatrix(8, 6)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,2);
	miniMatrix(8, 7)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(1,2); // there is a problem here
	miniMatrix(8, 8)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.numberOfBeads
			   - dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_product_sum(0,0)+hinge1.r_product_sum(1,1));
	miniMatrix(8,9)= +hinge1.r(1);
	miniMatrix(8,10)= -hinge1.r(0);

	// this is correct
	//cout << "minmatrix (9,8) " << miniMatrix(8, 7)<< endl;
	/*
	cout<< " minimatrix (9,4) " << miniMatrix(8,3) << " (9,5) " <<miniMatrix(8,4);
	cout<< " (9,7) " << miniMatrix(8,6 ) << " (9,8)" << miniMatrix(8,7);
	cout<< " (9,9) " << miniMatrix(8,8 ) << " (9,10)" << miniMatrix(8,9);
	cout<< " (9,11) " << miniMatrix(8,10)<<endl;
	*/
	///////////////////////////////

	miniVector(8)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(2)
           -+dragCoefficientVelocity * hinge1.averageviscosity * (hinge1.r_times_u_sum(0,1)-hinge1.r_times_u_sum(1,0))
           -hinge1.torque(2)- hinge1.excluded_volume_torque(0);// add excluded volume torque 
	// correct result
	//cout<< " Minivector (9) " << miniVector(8)<<endl;

	if (timeStep % TIME_STEP_DEBUG == 0)
	{
		cout<< " Hinge  Force" << hinge1.exluded_volume_force.transpose() <<endl;
		cout<< " Hinge  torque" << hinge1.excluded_volume_torque.transpose() <<endl;
	}

}



void miniMatrixAssembly(MatrixXd & miniMatrix, VectorXd &miniVector, const Hinge & hinge1, Hinge hinge2, double dragCoefficientVelocity, double dragCoefficientOmega){


	// when comparing to fortran code remember that indexing is 0 based here as oposed to fortran
	miniMatrix.fill(0);
	miniVector.fill(0);

	////////////////////
	miniMatrix(0,3) = 1;
	miniMatrix(0,7) = hinge1.r(2); // this is in Z
	miniMatrix(0,8) = -hinge1.r(1); // this is in Y
	miniMatrix(0,12) = -1;
	// I checked the numbers
	

	////////////////////
	miniMatrix(1,4) = 1;
	miniMatrix(1,6) = -hinge1.r(2); // this is in Z
	miniMatrix(1,8) = hinge1.r(0); // this is in X
	miniMatrix(1,13) = -1;
	// I checked the numbers

	
	////////////////////
	miniMatrix(2,5) = 1;
	miniMatrix(2,6) = hinge1.r(1); // this is in Z
	miniMatrix(2,7) = -hinge1.r(0); // this is in X
	miniMatrix(2,14) = -1;
	// I checked the numbers

	

	///////////////////////
	miniMatrix(3, 0)=  1;
	miniMatrix(3, 3)= -dragCoefficientVelocity* hinge1.averageviscosity * hinge1.numberOfBeads  ;
	miniMatrix(3, 7)= -dragCoefficientVelocity* hinge1.averageviscosity* hinge1.r_sum(2);
	miniMatrix(3, 8)= +dragCoefficientVelocity *hinge1.averageviscosity* hinge1.r_sum(1);
	miniMatrix(3,9)= -1;

	// I checked the numbers


	/******************/
	miniVector(3)= -dragCoefficientVelocity * hinge1.averageviscosity *  hinge1.fluid_velocity_sum(0) - hinge1.exluded_volume_force(0); // add excluded volume force
	// the number it returns is -0, I dont know if there might be any numerical implication to this. 

	

	/*
	cout<< " minimatrix (4,1) " << miniMatrix(3,0) << " (4,4) " <<miniMatrix(3,3);
	cout<< " (4,8) " << miniMatrix(3,7 ) << " (4,9)" << miniMatrix(3,8);
	cout<< " (4,10) " << miniMatrix(3,9 )<<endl;
	*/
	///////////////////

	miniMatrix(4, 1)=  1;
	miniMatrix(4, 4)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.numberOfBeads;
	miniMatrix(4, 6)= +dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(2);
	miniMatrix(4, 8)= -dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(0);
	miniMatrix(4,10)= -1;
	// correct result
	/*
	cout<< " minimatrix (5,2) " << miniMatrix(4,1) << " (5,5) " <<miniMatrix(4,4);
	cout<< " (5,7) " << miniMatrix(4,6 ) << " (5,9)" << miniMatrix(4,8);
	cout<< " (5,11) " << miniMatrix(4,10 )<<endl;
	*/
	//////////////////

	miniVector(4)= -dragCoefficientVelocity * hinge1.averageviscosity *  hinge1.fluid_velocity_sum(1)- hinge1.exluded_volume_force(1); // add excluded volume force
	// result is correct
	///////////////////

	miniMatrix(5, 2)=  1;
	miniMatrix(5, 5)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.numberOfBeads;
	miniMatrix(5, 6)= -dragCoefficientVelocity*hinge1.averageviscosity* hinge1.r_sum(1);
	miniMatrix(5, 7)= +dragCoefficientVelocity*hinge1.averageviscosity*hinge1.r_sum(0);
	miniMatrix(5,11)= -1;
	// result is correct
	/*
	cout<< " minimatrix (6,3) " << miniMatrix(5,2) << " (6,6) " <<miniMatrix(5,5);
	cout<< " (6,7) " << miniMatrix(5,6 ) << " (6,8)" << miniMatrix(5,7);
	cout<< " (6,12) " << miniMatrix(5,11 )<<endl;
	*/
	///////////////////////

	miniVector(5)= -dragCoefficientVelocity* hinge1.averageviscosity *  hinge1.fluid_velocity_sum(2)- hinge1.exluded_volume_force(2); // add excluded volume force

#ifdef PRINT_DEBUG_INFO
	cout<< "minivector(5) "<<miniVector(5)<<endl;
	cout<< "dragCoefficientVelocity "<<dragCoefficientVelocity<<endl;
	cout<< "hinge1.averageviscosity "<<hinge1.averageviscosity<<endl;
	cout<< "hinge1.fluid_velocity_sum(2) "<<hinge1.fluid_velocity_sum(2)<<endl;
	cout<< "hinge1.excludedVolumeForce(2) "<<hinge1.exluded_volume_force(2)<<endl;
	cout<< "hinge1.excludedVolumeForce "<<hinge1.exluded_volume_force<<endl;
#endif
	// again -0 as the answer
	
	///////////////////////

	miniMatrix(6,6)=1; // again!!!??? (6,6) is repeated
	miniMatrix(6, 4)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(2);
	miniMatrix(6, 5)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(1);
	miniMatrix(6, 6)= -dragCoefficientOmega * hinge1.averageviscosity *  hinge1.numberOfBeads  
					- dragCoefficientVelocity *hinge1.averageviscosity  *  ( hinge1.r_product_sum(2,2)+hinge1.r_product_sum(1,1)) ;           
	miniMatrix(6, 7)=  dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,1);
	miniMatrix(6, 8)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,2);
	miniMatrix(6,10)= +hinge1.r(2);
	miniMatrix(6,11)= -hinge1.r(1);
	// results are correct
	/*
	cout<< " minimatrix (7,7) " << miniMatrix(6,6) << " (7,5) " <<miniMatrix(6,4);
	cout<< " (7,6) " << miniMatrix(6,5 ) << " (7,7)" << miniMatrix(6,6);
	cout<< " (7,8) " << miniMatrix(6,7 ) << " (7,9)" << miniMatrix(6,8);
	cout<< " (7,11) " << miniMatrix(6,10 ) << " (7,12)" << miniMatrix(6,11)<<endl;
	*/
	/////////////////////////
	
	miniVector(6)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(0) 
		+dragCoefficientVelocity * hinge1.averageviscosity  * ( hinge1.r_times_u_sum(2,1)- hinge1.r_times_u_sum(1,2))
		-hinge1.torque(0) - hinge1.excluded_volume_torque(0);  // add excluded volume force 
	// correct result

	//////////////////////////

	miniMatrix(7, 3)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(2);
	miniMatrix(7, 5)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(0);
	miniMatrix(7, 6)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,1);
	miniMatrix(7, 7)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.numberOfBeads
           -dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_product_sum(2,2)+hinge1.r_product_sum(0,0));
	miniMatrix(7, 8)=  dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(1,2);
	miniMatrix(7,9)= -hinge1.r(2);
	miniMatrix(7,11)= +hinge1.r(0);
	// correct result 
	/*
	cout<< " minimatrix (8,4) " << miniMatrix(7,3) << " (8,6) " <<miniMatrix(7,5);
	cout<< " (8,7) " << miniMatrix(7,6 ) << " (8,8)" << miniMatrix(7,7);
	cout<< " (8,9) " << miniMatrix(7,8 ) << " (8,10)" << miniMatrix(7,9);
	cout<< " (8,12) " << miniMatrix(7,11)<<endl;
	*/
	///////////////////////////

	miniVector(7)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(1)
		+dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_times_u_sum(0,2)-hinge1.r_times_u_sum(2,0))
           -hinge1.torque(1)- hinge1.excluded_volume_torque(0);// add excluded volume torque 
	/*
	cout<<" dragCoefficientOmega " << dragCoefficientOmega;
	cout<<" averageviscosity " << hinge1.averageviscosity;
	cout<<" sumOmegaFluid(1) " << hinge1.sumOmegaFluid(1);
	cout<<" r_times_u_sum(0,2) " << hinge1.r_times_u_sum(0,2);
	cout<<" r_times_u_sum(2,0) " << hinge1.r_times_u_sum(2,0);
	cout<<" minivector(7) " << miniVector(7);
	cout<<endl;
	*/
	// correct result
	
	//////////////////////////////

	miniMatrix(8, 3)=  dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(1);
	miniMatrix(8, 4)= -dragCoefficientVelocity * hinge1.averageviscosity * hinge1.r_sum(0);
	miniMatrix(8, 6)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(0,2);
	miniMatrix(8, 7)= +dragCoefficientVelocity * hinge1.averageviscosity  *  hinge1.r_product_sum(1,2); // there is a problem here
	miniMatrix(8, 8)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.numberOfBeads
			   - dragCoefficientVelocity * hinge1.averageviscosity  *(hinge1.r_product_sum(0,0)+hinge1.r_product_sum(1,1));
	miniMatrix(8,9)= +hinge1.r(1);
	miniMatrix(8,10)= -hinge1.r(0);

	// this is correct
	//cout << "minmatrix (9,8) " << miniMatrix(8, 7)<< endl;
	/*
	cout<< " minimatrix (9,4) " << miniMatrix(8,3) << " (9,5) " <<miniMatrix(8,4);
	cout<< " (9,7) " << miniMatrix(8,6 ) << " (9,8)" << miniMatrix(8,7);
	cout<< " (9,9) " << miniMatrix(8,8 ) << " (9,10)" << miniMatrix(8,9);
	cout<< " (9,11) " << miniMatrix(8,10)<<endl;
	*/
	///////////////////////////////

	miniVector(8)= -dragCoefficientOmega * hinge1.averageviscosity * hinge1.sumOmegaFluid(2)
           -+dragCoefficientVelocity * hinge1.averageviscosity * (hinge1.r_times_u_sum(0,1)-hinge1.r_times_u_sum(1,0))
           -hinge1.torque(2)- hinge1.excluded_volume_torque(0);// add excluded volume torque 
	// correct result
	//cout<< " Minivector (9) " << miniVector(8)<<endl;

}



// Here is where the assembly and solution of the matrix for each fiber takes place, not really
// this is how daniel has it

void fiberMatrixAssembler(fiber & fiber1, double dragCoefficientVelocity, double dragCoefficientOmega){


	#ifdef PRINT_DEBUG_INFO
		

		for (int j = 0; j < fiber1.numberOfHinges; j++)
		{
			cout<< " Segment " << j << endl;
			cout<< " excluded volume force " << fiber1.hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fiber1.hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fiber1.hinges[j].excluded_volume_torque<<endl;
		}

	#endif




	// the system is A x = B
	// the size of A is the number of segments in the fiber times 9, as there are 9 equations per segment
	// number of segments is equal to the number of hinges minus one
	int xSize = 9*(fiber1.numberOfHinges-1);
	//MatrixXd A(xSize,xSize);
	MatrixXd A = MatrixXd::Zero(xSize,xSize);
	//VectorXd X(xSize);
	VectorXd X = VectorXd::Zero(xSize);
	//VectorXd B(xSize);
	VectorXd B= VectorXd::Zero(xSize);

	//cout<< " X " << X <<endl;

	//TO DO initialize this better
	//A.fill(0);
	//X.fill(0);
	//B.fill(0);
	//cout<< " matrix size "<< xSize << endl; 

	// miniMatrix and miniVector are the fraction of the matrix and vector related to that especific segment, 
	// this are computed and then are added to the biger matrix A
	MatrixXd miniMatrix(9,15);
	VectorXd miniVector(9);
	

	if (fiber1.numberOfHinges !=2 )
	{
	
		for (int i = 0; i < fiber1.numberOfHinges-1; i++)
		{
			miniMatrixAssembly(miniMatrix,miniVector,fiber1.hinges[i], fiber1.hinges[i+1],dragCoefficientVelocity,dragCoefficientOmega);
			if (i ==0)
			{
				copyMatrixRange(A,0,8,0,11,miniMatrix,0,3); // equivalent to Amat(1:9,1:12)=mat(1:9,4:15) in Daniel's code

			}else if (i== fiber1.numberOfHinges-2)
			{
				// first line
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)-3 << " upper bound column " << 9*(i)+5 << endl;

				//second line
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)+6 << " upper bound column " << 9*(i)+8 << endl;

				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)-3 , 9*(i)+5 , miniMatrix , 0 , 0);
				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)+6 , 9*(i)+8 , miniMatrix , 0 , 12);

			}else
			{
				//Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+12)= mat
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)-3 << " upper bound column " << 9*(i)+11 << endl;
				
				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)-3 , 9*(i)+11 , miniMatrix , 0 , 0);
				
			}

		// bvec(9*(i-1)+1: 9*i,1)=vec(1:9,1)
			copyVectorRange(B, 9*i, 9*(i+1)-1 , miniVector , 0); 

		}
	}else{
		// fiber has just one segment 
		miniMatrixAssembly(miniMatrix,miniVector,fiber1.hinges[0], fiber1.hinges[1],dragCoefficientVelocity,dragCoefficientOmega);
		copyMatrixRange(A,0,8,0,5,miniMatrix,0,3); // equivalent to Amat(1:9,1:6)=mat(1:9,4:9) in Daniel's code
		copyMatrixRange(A,0,8,6,8,miniMatrix,0,12); // equivalent to Amat(1:9,7:9) =mat(1:9,13:15) in Daniel's code
		copyVectorRange(B, 0, 8 , miniVector , 0); // equivalent to bvec(:,1)=vec(:,1)

		
	}


	//cout<<"minimatrix "<< miniMatrix<<endl;


	//cout << "time spent assemblying matrix  " << getMilliSpan(time1)<<endl;
	#ifdef PRINT_DEBUG_INFO
	cout<<"dragCoefficientVelocity " << dragCoefficientVelocity<< endl;
	cout<<"dragCoefficientOmega " << dragCoefficientOmega<< endl;
	cout<< " A " <<endl << A <<endl;
	cout<< " B " <<endl << B <<endl;
	#endif	
	X = A.colPivHouseholderQr().solve(B); // this vector has the velocities 
	#ifdef PRINT_DEBUG_INFO	
	cout<< " X " <<endl << X <<endl;
	//cout <<  " velocity result " << X; 
	#endif	


	if (fiber1.numberOfSegments == 1)
	{

		fiber1.hinges[0].velocity(0) = X(0);
		fiber1.hinges[0].velocity(1) = X(1);
		fiber1.hinges[0].velocity(2) = X(2);

		fiber1.hinges[1].velocity(0) = X(6);
		fiber1.hinges[1].velocity(1) = X(7);
		fiber1.hinges[1].velocity(2) = X(8);

	}else{

		for (int i = 0; i < fiber1.numberOfSegments; i++)
	{
		// X stores 9 degrees of freedom per hinge, the first 3 are the velocities in x, y and z
		// for the las hinge its values are stored in the last 3 positions of X
		fiber1.hinges[i].velocity(0) = X(9*i);
		fiber1.hinges[i].velocity(1) = X(9*i+1);
		fiber1.hinges[i].velocity(2) = X(9*i+2);
	}

	fiber1.hinges[fiber1.numberOfHinges-1].velocity(0)=  X(9*(fiber1.numberOfHinges-1) -3);
	fiber1.hinges[fiber1.numberOfHinges-1].velocity(1)=  X(9*(fiber1.numberOfHinges-1) -2);
	fiber1.hinges[fiber1.numberOfHinges-1].velocity(2)=  X(9*(fiber1.numberOfHinges-1) -1);
	
	}

}


//for debugging purposes

void fiberMatrixAssembler(fiber & fiber1, double dragCoefficientVelocity, double dragCoefficientOmega, int timeStep){


	#ifdef PRINT_DEBUG_INFO
		

		for (int j = 0; j < fiber1.numberOfHinges; j++)
		{
			cout<< " Segment " << j << endl;
			cout<< " excluded volume force " << fiber1.hinges[j].exluded_volume_force<<endl;
			cout<< " torque " << fiber1.hinges[j].torque<<endl;
			cout<< " excluded volume torque " << fiber1.hinges[j].excluded_volume_torque<<endl;
		}

	#endif




	// the system is A x = B
	// the size of A is the number of segments in the fiber times 9, as there are 9 equations per segment
	// number of segments is equal to the number of hinges minus one
	int xSize = 9*(fiber1.numberOfHinges-1);
	//MatrixXd A(xSize,xSize);
	MatrixXd A = MatrixXd::Zero(xSize,xSize);
	//VectorXd X(xSize);
	VectorXd X = VectorXd::Zero(xSize);
	//VectorXd B(xSize);
	VectorXd B= VectorXd::Zero(xSize);

	//cout<< " X " << X <<endl;

	//TO DO initialize this better
	//A.fill(0);
	//X.fill(0);
	//B.fill(0);
	//cout<< " matrix size "<< xSize << endl; 

	// miniMatrix and miniVector are the fraction of the matrix and vector related to that especific segment, 
	// this are computed and then are added to the biger matrix A
	MatrixXd miniMatrix(9,15);
	VectorXd miniVector(9);
	

	if (fiber1.numberOfHinges !=2 )
	{
	
		for (int i = 0; i < fiber1.numberOfHinges-1; i++)
		{
			miniMatrixAssembly(miniMatrix,miniVector,fiber1.hinges[i], fiber1.hinges[i+1],dragCoefficientVelocity,dragCoefficientOmega);
			if (i ==0)
			{
				copyMatrixRange(A,0,8,0,11,miniMatrix,0,3); // equivalent to Amat(1:9,1:12)=mat(1:9,4:15) in Daniel's code

			}else if (i== fiber1.numberOfHinges-2)
			{
				// first line
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)-3 << " upper bound column " << 9*(i)+5 << endl;

				//second line
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)+6 << " upper bound column " << 9*(i)+8 << endl;

				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)-3 , 9*(i)+5 , miniMatrix , 0 , 0);
				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)+6 , 9*(i)+8 , miniMatrix , 0 , 12);

			}else
			{
				//Amat(9*(i-1)+1: 9*i, 9*(i-1)-2: 9*(i-1)+12)= mat
				//cout<< " lower bound row " << 9*(i) << " upper bound row " << 9*(i+1)-1 << endl;
				//cout<< " lower bound column " <<  9*(i)-3 << " upper bound column " << 9*(i)+11 << endl;
				
				copyMatrixRange(A, 9*(i) , 9*(i+1)-1 , 9*(i)-3 , 9*(i)+11 , miniMatrix , 0 , 0);
				
			}

		// bvec(9*(i-1)+1: 9*i,1)=vec(1:9,1)
			copyVectorRange(B, 9*i, 9*(i+1)-1 , miniVector , 0); 

		}
	}else{
		// fiber has just one segment 
		miniMatrixAssembly(miniMatrix,miniVector,fiber1.hinges[0], fiber1.hinges[1],dragCoefficientVelocity,dragCoefficientOmega);
		copyMatrixRange(A,0,8,0,5,miniMatrix,0,3); // equivalent to Amat(1:9,1:6)=mat(1:9,4:9) in Daniel's code
		copyMatrixRange(A,0,8,6,8,miniMatrix,0,12); // equivalent to Amat(1:9,7:9) =mat(1:9,13:15) in Daniel's code
		copyVectorRange(B, 0, 8 , miniVector , 0); // equivalent to bvec(:,1)=vec(:,1)

		
	}


	//cout<<"minimatrix "<< miniMatrix<<endl;


	//cout << "time spent assemblying matrix  " << getMilliSpan(time1)<<endl;
	#ifdef PRINT_DEBUG_INFO
	cout<<"dragCoefficientVelocity " << dragCoefficientVelocity<< endl;
	cout<<"dragCoefficientOmega " << dragCoefficientOmega<< endl;
	cout<< " A " <<endl << A <<endl;
	cout<< " B " <<endl << B <<endl;
	#endif	
	X = A.colPivHouseholderQr().solve(B); // this vector has the velocities 
	#ifdef PRINT_DEBUG_INFO	
	cout<< " X " <<endl << X <<endl;
	//cout <<  " velocity result " << X; 
	#endif	


	if (fiber1.numberOfSegments == 1)
	{

		fiber1.hinges[0].velocity(0) = X(0);
		fiber1.hinges[0].velocity(1) = X(1);
		fiber1.hinges[0].velocity(2) = X(2);

		fiber1.hinges[1].velocity(0) = X(6);
		fiber1.hinges[1].velocity(1) = X(7);
		fiber1.hinges[1].velocity(2) = X(8);

	}else{

		for (int i = 0; i < fiber1.numberOfSegments; i++)
	{
		// X stores 9 degrees of freedom per hinge, the first 3 are the velocities in x, y and z
		// for the las hinge its values are stored in the last 3 positions of X
		fiber1.hinges[i].velocity(0) = X(9*i);
		fiber1.hinges[i].velocity(1) = X(9*i+1);
		fiber1.hinges[i].velocity(2) = X(9*i+2);
	}

	fiber1.hinges[fiber1.numberOfHinges-1].velocity(0)=  X(9*(fiber1.numberOfHinges-1) -3);
	fiber1.hinges[fiber1.numberOfHinges-1].velocity(1)=  X(9*(fiber1.numberOfHinges-1) -2);
	fiber1.hinges[fiber1.numberOfHinges-1].velocity(2)=  X(9*(fiber1.numberOfHinges-1) -1);
	
	}

}



//The fiber motion matrix is assembled and the system of equations is solved
void solveFiberMotion( vector<fiber> &fibers, double fiberRadius){

	double dragCoefficientVelocity = 6 * pi * fiberRadius;
	double dragCoefficientOmega = 8 * pi * pow(fiberRadius,3);



	//#pragma omp parallel for
	for (int i = 0; i < fibers.size(); i++)
	{
		fiberMatrixAssembler(fibers[i],dragCoefficientVelocity,dragCoefficientOmega);

	}



}

// for debugging purposes
void solveFiberMotion( vector<fiber> &fibers, double fiberRadius, int timeStep){

	double dragCoefficientVelocity = 6 * pi * fiberRadius;
	double dragCoefficientOmega = 8 * pi * pow(fiberRadius,3);

	


	
	
	#pragma omp parallel for
	for (int i = 0; i < fibers.size(); i++)
	{
		if (!fibers[i].hinges[0].isStationary)
		{

		
		//cout<<" number of maximum threads "<<omp_get_max_threads()<<", threads being  used :" << omp_get_num_threads()<<endl;
		fiberMatrixAssembler(fibers[i],dragCoefficientVelocity,dragCoefficientOmega);
		}
	}



}