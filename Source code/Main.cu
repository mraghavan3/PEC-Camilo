// This version works with comsol files
#pragma once
#include <map>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream> 

#define _USE_MATH_DEFINES // Math constants
#include <math.h>
#include "variables.h"
#include "bending torque.h"
#include "check fiber breakage.h"
#include "SEgments parameters computation.h"
#include "Fibers Motion Solver.h"
#include "utilities.h"
#include "Find Neighbors.cuh"
#include "Excluded Volume.h"
#include "InputData.h"
#include "time integration.h"
#include "velocityMeshUtilities.h"
#define PRINT_DEBUG_INFO_GLOBAL_TIMING
//#define PRINT_DEBUG_INFO
#define IS_WINDOWS 0
#define IS_LINUX 1
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
const int sytemOS = IS_WINDOWS;
#else
const int sytemOS = IS_LINUX;
#endif


// This is a direct translation of the mechanistic Model from Fortran to C++
// Started june/2013
using namespace std;
using namespace Eigen;


	
void initializeNamelist(){
	s_mapVariablesValues["r_bead" ]=fiber_radius_name;
	s_mapVariablesValues["viscosity" ]=viscosity_name;
	s_mapVariablesValues["dt" ]=dt_name;
	s_mapVariablesValues["nbr_intgr" ]=number_integrations_name;
	s_mapVariablesValues["writ_freq" ]=writing_frequency_name;
	s_mapVariablesValues["E_Young" ]=E_young_name;
	s_mapVariablesValues["L_seg" ]=L_seg_name;
	s_mapVariablesValues["max_alpha" ]=critical_angle_name;
	s_mapVariablesValues["break_freq" ]=check_breakage_frequency_name;
	s_mapVariablesValues["gamma_dot" ]=gamma_dot_name;
}
int readData(){

	// Approach to read namelists as fortran does but way more primitive

	// need a list of values that represent the variables names
	//string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Fibers.in"); // for folders use "/" do not use"\"
	
	  if (tempFile.is_open())
	  {
		  cout<<" Reading input parameters \n";
		  initializeNamelist();
		  while(getline (tempFile,line))  //read line 
		  {
			  size_t equalPosition = line.find('=');
			  if (equalPosition< line.length()) // gets rid of lines that do not have equals
			  {
				  line.replace(equalPosition, 1 , " "); //delete equal so no more steps are needed after tokenization
			  stringstream lineTokens(line); //Tokenization of line
			  string variableName;
			  lineTokens>>variableName; // first token is the variable name
			  cout<<" "<< variableName << " ";
			  switch (s_mapVariablesValues[variableName]) // look for the variable name "value" in our map and switch it
			  {
				  case gamma_dot_name:
					  lineTokens >>gamma_dot ; 
					cout<< " gamma_dot is : " << gamma_dot<< "  ";
				  break;

				   case check_breakage_frequency_name:
					  lineTokens >>check_breakage_frequency ; 
					cout<< " check_breakage_frequency is : " << check_breakage_frequency<< "  ";
				  break;

					case critical_angle_name:
					  lineTokens >>critical_angle ; 
					cout<< " critical_angle is : " << critical_angle<< "  ";
				  break;
				  case L_seg_name:
					  lineTokens >>L_segment ; 
					cout<< " L_seg is : " << L_segment<< "  ";
				  break;

				  case writing_frequency_name:
					lineTokens >>writing_frequency ; 
					cout<< " writing_frequency is : " << writing_frequency<< "  ";
				  break;

				   case number_integrations_name:
					lineTokens >>number_integrations ; 
					cout<< " number_integrations is : " << number_integrations << "  ";
				  break;

				  case dt_name:
					lineTokens >>dt ; 
					cout<< " dt is : " << dt << "  ";
				  break;

				  case viscosity_name:
					lineTokens >>viscosity ; 
					cout<< " viscosity is : " << viscosity << "  ";
				  break;

				case fiber_radius_name:
					lineTokens >>fiber_radius ; 
					cout<< " fiber_radius is : " << fiber_radius << "  ";
				  break;
			  
			  case E_young_name:
					lineTokens >>E_young ; 
					cout<< " Eyoung is : " << E_young << "  ";
				  break;
			  
			  default:
				  break;
			  }
			
			  cout<<endl;

			  }
			
		  }

		  tempFile.close();
	 }

	return 0;
}


int main(){
	

	// Initialization of the program, reading input files
	
	vector<fiber> fibers;
	
	cout<< " Mechanistic model of fibers moving in a fluid \n \n" ; 
	
	// Variables definition
	ofstream outputFile;
	ofstream framesFile;
	thrust::host_vector<int> fiberIndices;
	thrust::host_vector<int> segmentIndices;
	thrust::host_vector<long long> possibleCollisions;

	// Variables needed for 3D transient mesh flow field
	vector<double> velMeshPoints;    //Stores coordinates of all vertices in the mesh as follows
									 //X0,X1,...Xn, Y0,Y1,...,Yn, Z0,Z1,..Zn

	vector<int> velMeshConnectivity; // Stores the vertices that make up the tetrahedra
									 // each component is an node index, there are 4 index per tet
									 // N00,N01,N02,N03,N10,N11,N12,N13,...,Nn0,Nn1,n2,Nn3
	
	vector<int> velMeshNeighbors;	 // Stores the face-neighbors for each tetrahedron
									 // each component is a tet index, there are 4 index per tet
									 // N00,N01,N02,N03,N10,N11,N12,N13,...,Nn0,Nn1,n2,Nn3
									 // 3 for Tri , 6 for hex and so on
	
	vector<int> tetraHingeID;		// Stores the tet index where a certain hinge is
										
	vector<double> velTime1;		// Stores the velocity of all the nodes at a certaint time
									//Vx0,Vx1,...Vxn, Vy0,Vy1,...,Vyn, Vz0,Vz1,..Vzn
	
	vector<double> velTime2;		//same as velTime1

	vector<double> meshVorticity;		// value of vorticity at different nodes

	vector<int> meshCorners;        // in 2D nodesID of nodes that are corners, We will add statci fibers along those points

	vector<int> boundaryEdgeIndex;   // for each element stores -1 if there are no boundary edges and if there is one stores which face is boundary



	// Read Data
	outputFile.open("Output/data2.in",ios::out);
	
	readBorder(border_coords, "Input/concise_marco.txt"); // wall information is stored in border coords. TODO should be able to receive 2D or 3D
	readFlowFieldData(flowField, "Input/Coords.txt");     // Velocity field info, stored in flow field. TODO 
	readData();											  // information about the fibers 
	
	readInitialFibersPositions(fibers, hinges);// self explanatory 
	
	int nodesPerElement;

	readMeshFileTri("Input/rib_mesh.mphtxt",velMeshConnectivity,velMeshPoints,meshCorners,nodesPerElement);
	readVelFile2D("Input/fields_in_rib2.txt", velTime1,meshVorticity);
	int numNodes = velTime1.size()/3;
	int numElements = velMeshConnectivity.size()/nodesPerElement; // 6 for triangles

	generateCornerFibers(meshCorners,fibers,velMeshPoints,0.1,numNodes);
	cout<<endl<<endl<<" Simulation summary:  "<<endl;
	cout<<" Number of fibers : " << fibers.size() << " \n";
	cout<<" Number of elements : " << numElements << " \n";
	cout<<" Number of numNodes : " << numNodes << " \n";
	cout<<" Number of nodes per element : " << nodesPerElement << " \n";
	// Process the mesh
	
	velMeshNeighbors.resize(numElements*3);//3 for triangles TODO get this value from mesh Type
	

	

	boundaryEdgeIndex.resize(numElements);//assuming there is just 1 boundary per element
	find_tri_neighbors(numElements,nodesPerElement,velMeshConnectivity,velMeshNeighbors); // find neighboring elements for each element
	int elementTotest = 19;
	cout<< " the neighbors of element " <<elementTotest<< " are : "<<velMeshNeighbors[elementTotest*3+0 ]<< " "<<velMeshNeighbors[elementTotest*3+1 ]<< " "<<velMeshNeighbors[elementTotest*3+2 ]<<endl;
	
	findBoundaryEdges(numElements,3,velMeshNeighbors,boundaryEdgeIndex); //finds boundary edges indices for each element, helpful for excluded volume forces

	find_Initial_tetID(velMeshPoints,velMeshConnectivity,fibers,numElements,numNodes,2,nodesPerElement);//meshtype 2, tri elements
	

	/* initialize fiber and segment indices;
	*/
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


	double inertia_moment = (pi/4)*pow(fiber_radius,4);
	double * boundingBoxes;
	boundingBoxes = (double*) malloc(sizeof(double)* hinges.size() *6 );

	//main loop goes here, right now it is just one integration
	int timeA =getMilliCount();
	int timeB;

	int frames = 0;

	#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
	number_integrations=5;
	#endif
	#ifdef PRINT_DEBUG_INFO
	number_integrations=1;
	#endif

	//number_integrations =1; // Delete this line when not debugging
	for (int step = 0; step < number_integrations; step++)
	{
		if (step% check_breakage_frequency == 0 )
		{

			bendingTorqueForAllFibers(fibers,hinges,E_young, inertia_moment); // TO DO, Hinges can be removed this is bending_torque_whole
			checkFibersBreakage(fibers, critical_angle, fiberIndices, segmentIndices);         //this is  fiber_damage
			#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
			timeA =getMilliCount();
				#endif

			find_Initial_tetID(velMeshPoints,velMeshConnectivity,fibers,numElements,numNodes,2,nodesPerElement);//meshtype 2, tri elements
			#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "find_Initial_tetID done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
				
				#endif
			if (fibers.size()!= 1)
			{
				#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
				#endif

				findNeighbours(fibers,fiber_radius, possibleCollisions, fiberIndices); // To do 
								
				if (possibleCollisions.size() !=0)
				{
					refineNeighbors(fibers, fiber_radius, possibleCollisions, fiberIndices, segmentIndices);
				}
				
				#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "Find Neighbors done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
				cout<<" possible collision size" << possibleCollisions.size()<<endl;
				#endif

			}
			
		}

		

		// 
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif
				 checkHingesElements(velMeshPoints,velMeshConnectivity,velMeshNeighbors,fibers,numElements,numNodes,2,nodesPerElement);//2 for tri

		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "checkHingesElements done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif


///////////////**************************//////////////			
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif

				 //fibers_parameters_computation(fibers,fiber_radius); // fiber_par_calc
				 fibers_parameters_computation_mesh(fibers,fiber_radius,velMeshPoints,velMeshConnectivity,velTime1,meshVorticity,numElements,numNodes,nodesPerElement); // 6 is quadratic triangles
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "fibers_parameters_computation done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif

///////////////**************************//////////////			
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif

		//calculateFiberInteractions(fibers,fiber_radius,possibleCollisions,fiberIndices,segmentIndices,step); // prints debugging info
		calculateFiberInteractions(fibers,fiber_radius,possibleCollisions,fiberIndices,segmentIndices);

		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "calculateFiberInteractions done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif

#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
#endif
				 interactionsWallsTri(fibers,numNodes,velMeshPoints,velMeshConnectivity,velMeshNeighbors,fiber_radius,100,10*fiber_radius,nodesPerElement);
   #ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "interactionsWallsTri done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif

///////////////**************************//////////////	
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif

		bendingTorqueForAllFibers(fibers, hinges, E_young, inertia_moment);
		
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "bendingTorqueForAllFibers done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif

///////////////**************************//////////////	
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif

			solveFiberMotion(fibers,fiber_radius,step); // mot

		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "solveFiberMotion done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif
		
///////////////**************************//////////////	
		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				 timeA =getMilliCount();
		#endif

		updatePositions(fibers,dt);
		

		#ifdef PRINT_DEBUG_INFO_GLOBAL_TIMING
				cout << "updatePositions done, time elpased:  " << getMilliSpan(timeA)<<" ms"<<endl;
	    #endif
///////////////**************************//////////////	

		#ifdef PRINT_DEBUG_INFO
		for (int j= 0; j < fibers[0].numberOfHinges; j++)
		{
			cout<<endl<<endl;
			cout<< " Hinge " << j<< " r_times_u_sum "<< fibers[0].hinges[j].r_times_u_sum<<endl;
			cout<< " Hinge " << j<< " fluid_velocity_sum "<< fibers[0].hinges[j].fluid_velocity_sum<<endl;
			cout<< " Hinge " << j<< " sumOmegaFluid "<< fibers[0].hinges[j].sumOmegaFluid<<endl;
			cout<< " Hinge " << j<< " averageviscosity "<< fibers[0].hinges[j].averageviscosity<<endl;
			cout<< " Hinge " << j<< " r_sum "<< fibers[0].hinges[j].r_sum<<endl;
			cout<< " Hinge " << j<< " r_product_sum "<< fibers[0].hinges[j].r_product_sum<<endl;
			cout<< " Hinge " << j<< " exluded_volume_force "<< fibers[0].hinges[j].exluded_volume_force<<endl;
			cout<< " Hinge " << j<< " excluded_volume_torque "<< fibers[0].hinges[j].excluded_volume_torque<<endl;
			cout<< " Hinge " << j<< " torque "<< fibers[0].hinges[j].torque<<endl;
			cout<< " Hinge " << j<< " velocity "<< fibers[0].hinges[j].velocity<<endl;
			cout<< " Hinge " << j<< " position "<< fibers[0].hinges[j].position<<endl;
			int elementID = fibers[0].hinges[j].elementID;
			cout<< " Hinge " << j << " Element " << elementID<<endl;
			cout<< " nodes "<< " : " ;
				for (int kk = 0; kk < nodesPerElement; kk++)
				{
					cout << velMeshConnectivity[elementID*nodesPerElement+kk]<<" ";
				}
				cout<<endl;
			



		}
		
		cout<<endl;
		#endif
		
		//cout<<" Time needed for updatePositions is : " << getMilliSpan(timeB) << "ms" << endl;
		//cout<<" Integration number:  " << step<<endl; 
		// Write output file
		if (step%writing_frequency == 0)
		{
			
			frames ++;
			printPositions(fibers,outputFile,frames);
			cout<<" frame number : " << frames<< " time elapsed :" <<getMilliSpan(timeA)/1000<< " s " <<" number of collisions : " << possibleCollisions.size()<<endl;
			cout<<" integration number " << step << " of " << number_integrations <<endl;
			/*for (int i = 0; i < possibleCollisions.size(); i++)
			{
					int fiberIndexA = fiberIndices[((long)(possibleCollisions[i]>>32)) ];
					int segmentA = ((long)(possibleCollisions[i]>>32)) - segmentIndices[fiberIndexA];
					int fiberIndexB = fiberIndices[(long) possibleCollisions[i]];
					int segmentB = (long) possibleCollisions[i] - segmentIndices[fiberIndexB];

			cout<<" fiber A " << fiberIndexA << " segment " << segmentA <<" Fiber B " <<fiberIndexB;
				cout<<" segment " << segmentB<<endl;

			}
			*/


			timeA= getMilliCount();
			/*
			Hinge HingeA0 = fibers[0].hinges[0];
			Hinge HingeA1 = fibers[0].hinges[1];
			Hinge HingeB0 = fibers[1].hinges[0];
			Hinge HingeB1 = fibers[1].hinges[1];
			
			cout<<" segment A  : " << HingeA0.position.transpose() << " - " << HingeA1.position.transpose()<<endl;
		 cout<<" segment B  : " << HingeB0.position.transpose() << " - " << HingeB1.position.transpose()<<endl;
		cout<<" segment A Excluded volume force : "<< HingeA0.exluded_volume_force.transpose() << endl;
		cout<<" segment B Excluded volume force : "<< HingeB0.exluded_volume_force.transpose() << endl;
		cout<<" segment A Excluded volume torque : "<< HingeA0.excluded_volume_torque.transpose() << endl;
		cout<<" segment B Excluded volume torque : "<< HingeB0.excluded_volume_torque.transpose() << endl;
		
		cout<<endl;*/
		}
		
	}
	

	free(boundingBoxes);
	outputFile.close();

	std::cin.get();


  return 0;
    

}
