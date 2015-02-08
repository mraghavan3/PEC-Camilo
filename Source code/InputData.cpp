#include "InputData.h"





//Initializes border_coords and populates it
int readBorder(std::vector<double>  & border_coords, string filename){

///////////////// Reading boundary segments ////////////////
	//The domain is enclosed by a 2D frame, written in concise_marco.txt
	int number_boundary_segments  ;
	std::string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/concise_marco.txt"); // for folders use "/" do not use"\"


  if (tempFile.is_open())
  {
      
	  getline (tempFile,line);
	  number_boundary_segments = stoi(line); // reads number of segments
	  border_coords.resize(number_boundary_segments*6,0); // "preallocates the size of the vector saves some assignments later

		for (int i = 0; i < number_boundary_segments; i++)
			{

	  		getline (tempFile,line);
			
			//cout<<endl;
			for (int j = 0; j < 6; j++)
				{
					if (j!=2 && j!=5 ) // marco right now is 2d so Z is 0, indices for z are i*6+2 and i*6+5
					{
						//border_coords.push_back(stod(line,&sz));
						border_coords[i*6+j] = stod(line,&sz);
						line = line.substr(sz);

					}
					//cout<<border_coords[i*6+j]<<" ";
				}

			 }

		tempFile.close();
		cout<<" Border data has been succesfully loaded"<<endl;
	}

  return 0;
}



	//Initializes flowField and populates it
int readVelocities(std::vector<double>  & flowField, string filename){
	string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Vels.txt"); // for folders use "/" do not use"\"


  if (tempFile.is_open())
  {
       getline (tempFile,line);
	  int number_of_mesh_nodes = stoi(line); // reads number of segments
	  
	  for (int i = 0; i < number_of_mesh_nodes; i++)
			{

	  		getline (tempFile,line);
			//cout<<stod(line,&sz)<< " ";
			//cout<<endl;
			for (int j = 3; j < 4; j++)// velocities indices are Vx = 3 Vy =4 Vz=5
				{
				//border_coords.push_back(stod(line,&sz));
				flowField[i*8+j] = stod(line,&sz);
				//double temp = stod(line,&sz);
				line = line.substr(sz);
				//cout << temp<<" ";
				
				}

			 }

		cout<< " Velocity field succesfully loaded"<<endl;
		tempFile.close();
	}

	return 0;
}


int readViscosity(std::vector<double>  & flowField, string filename){
		
	string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Viscosity.txt"); // for folders use "/" do not use"\"
	
  if (tempFile.is_open())
  {
      getline (tempFile,line);
	  int number_of_mesh_nodes = stoi(line); // reads number of nodes
	  
	  for (int i = 0; i < number_of_mesh_nodes; i++)
			{

	  		getline (tempFile,line);
			//cout<<stod(line,&sz)<< " ";
			//cout<<endl;
			int j = 7; //index of viscosity is 7
			//border_coords.push_back(stod(line,&sz));
			flowField[i*8+j] = stod(line,&sz);
			//double temp = stod(line,&sz);
			line = line.substr(sz);
			//cout << temp<<" ";
			//cout<<flowField[i*8+j]<<" ";
			
			 }

	  
	    cout<< " viscosity field succesfully loaded"<<endl;
		tempFile.close();
	 }
  	
  return 0;
}
int readVorticity(std::vector<double>  & flowField, string filename){
	string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Vorticity.txt"); // for folders use "/" do not use"\"


	  if (tempFile.is_open())
	  {
		  getline (tempFile,line);
		  int number_of_mesh_nodes = stoi(line); // reads number of nodes
	  
		  for (int i = 0; i < number_of_mesh_nodes; i++)
				{

	  			getline (tempFile,line);
				//cout<<stod(line,&sz)<< " ";
				//cout<<endl;
				int j = 6;  // Index of the vorticity is 6
				//border_coords.push_back(stod(line,&sz));
				flowField[i*8+j] = stod(line,&sz);
				//double temp = stod(line,&sz);
				line = line.substr(sz);
				//cout << temp<<" ";
				//cout<<flowField[i*8+j]<<" ";
				
				 }

		  cout<< " Vorticity field succesfully loaded"<<endl;
		  tempFile.close();
		  
		 }

		return 0;
	}

int readFlowFieldData(std::vector<double>  & flowField, string filename){
	// The Velocity, viscosity and vorticity fields can be found in 
	// coords.txt vels.txt vorticity.txt viscosity.txt
	// It would be nice if both 2d and 3d were supported. 

		// First flowField needs to be preallocated with the number of nodes 
		// and the positions x y and z are red as well
		
	string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Coords.txt"); // for folders use "/" do not use"\"


  if (tempFile.is_open())
  {
      
	  getline (tempFile,line);
	  int number_of_mesh_nodes = stoi(line); // reads number of nodes in mesh
	  flowField.resize(number_of_mesh_nodes*8,0); // "preallocates the size of the vector, saves having to make some assignments later

	  for (int i = 0; i < number_of_mesh_nodes; i++)
			{

	  		getline (tempFile,line);
			//cout<<stod(line,&sz)<< " ";
			//cout<<endl;
			
			for (int j = 0; j < 2; j++)
				{
					if (j!=2 && j!=5 ) // marco right now is 2d so Z is 0, indices for z are i*6+2 and i*6+5
					{
						//border_coords.push_back(stod(line,&sz));
						flowField[i*8+j] = stod(line,&sz);
						//double temp = stod(line,&sz);
						line = line.substr(sz);
						//cout << temp<<" ";

					}
					//cout<<flowField[i*8+j]<<" ";
				}

			 }

		tempFile.close();

		cout<< " Flow field node coordinates succesfully loaded" <<endl;
		readVelocities(flowField, "Input/Vels.txt");
		readVorticity(flowField, "Input/Vorticity.txt");
		readViscosity(flowField, "Input/Viscosity.txt");
	 }
  return 0;
 }

int readInitialFibersPositions(vector<fiber> &fibers, vector<Hinge> &hinges){
	string::size_type sz;     // alias of size_t
	string line;
	ifstream tempFile ("Input/Initial_Positions.txt"); // for folders use "/" do not use"\"
	fiber fiber1;
	Hinge hinge1;
	int hingesIndex=0;
	const int X=0;
	const int Y=1;
	const int Z=2;
  if (tempFile.is_open())
  {
      getline (tempFile,line);
	  int number_of_fibers = stoi(line); // reads number of fibers
	  
	 for (int i = 0; i < number_of_fibers; i++)
		//for (int i = 0; i < 1; i++)
	  {

	  		getline (tempFile,line);
			int number_of_hinges = stoi(line);
			fiber1.numberOfHinges = number_of_hinges;
			fiber1.numberOfSegments = number_of_hinges-1;
			fiber1.midpoints.resize(fiber1.numberOfSegments);
			fiber1.segments.resize(fiber1.numberOfSegments);
			fiber1.firstHinge = hingesIndex;
			fibers.push_back(fiber1);
			
			//cout<< "  fiber  number : "<< i << " number of hinges = " << fibers[i].numberOfHinges;
			//cout<<endl;
			

			for (int j = 0; j < number_of_hinges; j++)
			{
				hinge1.fiber_ID = i;
				getline (tempFile,line);
				hinge1.isStationary = stoi(line,&sz);
				line = line.substr(sz);
				hinge1.position(X)=stod(line,&sz);
				line = line.substr(sz);
				hinge1.position(Y)=stod(line,&sz);
				line = line.substr(sz);
				hinge1.position(Z)=stod(line,&sz);

				hinge1.sumOmegaFluid.fill(0);
				hinge1.fluid_velocity_sum.fill(0);
				hinge1.r_sum.fill(0);
				hinge1.r_product_sum.fill(0);
				hinge1.r_times_u_sum.fill(0);
				hinge1.exluded_volume_force.fill(0);
				hinge1.excluded_volume_torque.fill(0);
				hinge1.averageviscosity = 0;
				hinge1.velocity.fill(0);
				hinge1.vorticity.fill(0);
				hinge1.elementID =-1;

				hinges.push_back(hinge1 );

				fibers[i].hinges.push_back(hinge1);

				//cout<< " hinge : " << j << " coords :" << hinges[hingesIndex].position.x() << " " << hinges[hingesIndex].position.y() << " "<<hinges[hingesIndex].position.z();
				//cout<<endl;
				hingesIndex++;
			}
		 }
	    cout<< " initial positions loaded"<<endl;
		tempFile.close();
	 }
  cout<< " Number of hinges : "<< hingesIndex<<"\n";

  return 0;

}


void generateCornerFibers(vector<int > corners,vector<fiber>&fibers,vector<double> meshPoints, double depth, int numNodes  ){

	fiber fiber1;
	Hinge hinge1;
	int initialNumFibers = fibers.size();
	int number_of_hinges =2;
	int hingesIndex = fibers[fibers.size()].firstHinge + fibers[fibers.size()].numberOfHinges;
	
	for (int corner = 0; corner < corners.size(); corner++)
	{
		

		fiber1.numberOfHinges = number_of_hinges;
		fiber1.numberOfSegments = number_of_hinges-1;
		fiber1.midpoints.resize(fiber1.numberOfSegments);
		fiber1.segments.resize(fiber1.numberOfSegments);
		fiber1.firstHinge = hingesIndex;
		fibers.push_back(fiber1);

		//first hinge
		hinge1.fiber_ID = initialNumFibers+corner;
		hinge1.isStationary = true;
		hinge1.position(0)=meshPoints[ corners[corner] + numNodes * 0  ] ;
		hinge1.position(1)=meshPoints[ corners[corner] + numNodes * 1  ] ;
		hinge1.position(2)=meshPoints[ corners[corner] + numNodes * 2  ] ;
		hinge1.sumOmegaFluid.fill(0);
		hinge1.fluid_velocity_sum.fill(0);
		hinge1.r_sum.fill(0);
		hinge1.r_product_sum.fill(0);
		hinge1.r_times_u_sum.fill(0);
		hinge1.exluded_volume_force.fill(0);
		hinge1.excluded_volume_torque.fill(0);
		hinge1.averageviscosity = 0;
		hinge1.velocity.fill(0);
		hinge1.vorticity.fill(0);
		hinge1.elementID =-1;
		fibers[initialNumFibers+corner].hinges.push_back(hinge1);
		hingesIndex++;
		//second hinge

		hinge1.fiber_ID = initialNumFibers+corner;
		hinge1.isStationary = true;
		hinge1.position(0)=meshPoints[ corners[corner] + numNodes * 0  ] ;
		hinge1.position(1)=meshPoints[ corners[corner] + numNodes * 1  ] ;
		hinge1.position(2)=meshPoints[ corners[corner] + numNodes * 2  ] + depth;//add a distance = depth
		hinge1.sumOmegaFluid.fill(0);
		hinge1.fluid_velocity_sum.fill(0);
		hinge1.r_sum.fill(0);
		hinge1.r_product_sum.fill(0);
		hinge1.r_times_u_sum.fill(0);
		hinge1.exluded_volume_force.fill(0);
		hinge1.excluded_volume_torque.fill(0);
		hinge1.averageviscosity = 0;
		hinge1.velocity.fill(0);
		hinge1.vorticity.fill(0);
		hinge1.elementID =-1;
		fibers[initialNumFibers+corner].hinges.push_back(hinge1);
		hingesIndex++;


	}


}


int printPositions(vector<fiber> fibers, ofstream & outputFile, int frames) {

	if (outputFile.is_open())
	{
		outputFile << fibers.size()<<"\n";

		for (int i = 0; i < fibers.size(); i++)
		{
			outputFile << fibers[i].numberOfHinges<<"\n";
			for (int j = 0; j < fibers[i].numberOfHinges; j++)
			{
				outputFile << fibers[i].hinges[j].position.x() << " ";
				outputFile << fibers[i].hinges[j].position.y() << " ";
				outputFile << fibers[i].hinges[j].position.z() << " ";
				outputFile << "\n";
	     	}
		}
	}

	
	ofstream framesFile;
	framesFile.open("Output/nbr_frames.txt",ios::out);
	if (framesFile.is_open())
	{
		framesFile << frames;
	}
	framesFile.close();

return 0;
}

void readVelFile2D(string Filename, vector<double> &velValues, vector<double> &vorticity  ){

	// reads mesh info from a comsol file, 
	// reads coordinates of the nodes, and tthe connectivity of the elements, 
	// also has information about the corners

	stringstream ss;
	string line;
	string tempString;
	
	// openFile

	ifstream tempFile (Filename); // add includes for this

	if (tempFile.is_open())
	  {

		  // skip 4 lines
		  for (int i = 0; i < 4; i++)
		  {
			  getline (tempFile,line);
		  }

		  // read number of nodes in this line

		    getline (tempFile,line);
		  ss<<line;
		  int numNodes;
		  ss>>tempString; // %
		  ss>>tempString; // Nodes:
		  ss>>numNodes;
		  ss=stringstream();
		  velValues.resize(numNodes*3);
		  vorticity.resize(numNodes*3);

		  // skip 4 lines
		  for (int i = 0; i < 4; i++)
		  {
			  getline (tempFile,line);
		  }

		  // read values

		  for (int i = 0; i < numNodes; i++)
		  {
			  getline (tempFile,line);
		  ss<<line;
		  ss>> tempString; ss>>tempString; // first 2 values are x and y coords

		  ss>>velValues[i ]; // VelX
		  ss>>velValues[i+numNodes ]; // Vely

		  ss>>vorticity[i + numNodes*2 ]; // vorticity Z
		   ss=stringstream();
		  }

		  tempFile.close();
		  cout<<"velocity file read succesfull, nodes: "<<numNodes<<endl; 
	}
	else{

	cout<< " the file with vel info could not be opened" <<endl;
	}


}

void readMeshFileTri(string Filename, vector<int> &meshConnectivity, vector<double> &meshPoints, vector<int> & meshCorners, int & nodesPerElement   ){

	// reads mesh info from a comsol file, 
	// reads coordinates of the nodes, and tthe connectivity of the elements, 
	// also has information about the corners

	stringstream ss;
	string line;
	string tempSubstring;
	// openFile

	ifstream tempFile (Filename); // add includes for this

	if (tempFile.is_open())
	  {
		  // skip 19 lines
		  for (int i = 0; i < 18; i++)
		  {
			  getline (tempFile,line);
		  }

		  //read number of nodes

		  getline (tempFile,line);
		  ss<<line;
		  int numNodes;
		  ss>>numNodes;
		  ss=stringstream();
		  meshPoints.resize(numNodes*3);

		  
		  // skip 3 lines

		  for (int i = 0; i < 3; i++)
		  {
			  getline (tempFile,line);
		  }

		  // start reading mexh points coords

		  for (int i = 0; i < numNodes; i++)
		  {
			  getline (tempFile,line);
		  ss<<line;
		  
		  ss>>meshPoints[i ]; // xCoord
		  ss>>meshPoints[i +numNodes*1]; //ycoord
		  ss=stringstream();
		  
		  }

		  //skip 9 lines
		  for (int i = 0; i < 9; i++)
		  {
			  getline (tempFile,line);
			  
		  }

		  //read number of coorners

		  getline (tempFile,line);
		  
		  ss<<line;
		  
		  int numCorners;
		  ss>>numCorners;
		  ss=stringstream();
		  
		  meshCorners.resize(numCorners);

		  // skip a line
		  getline (tempFile,line);

		  //read corners IDs
		  for (int i = 0; i < numCorners; i++)
		  {
			  getline (tempFile,line);
		  ss<<line;
		  ss>>meshCorners[i ]; // xCoord
		   ss=stringstream();
		  }

		  // skip 5 lines

		  for (int i = 0; i < 5; i++)
		  {
			  getline (tempFile,line);
			  
		  }

		  // // read number of geometric entity indices
		  getline (tempFile,line);
		  ss<<line;
		  int linesToSkip;
		  ss>>linesToSkip;
		  ss=stringstream();

		  // skip that number + 11
		  for (int i = 0; i < linesToSkip+11; i++)
		  {
			  getline (tempFile,line);
		  }

		  // // read number of edges
		  getline (tempFile,line);
		  ss<<line;
		  linesToSkip;
		  ss>>linesToSkip;
		  ss=stringstream();

		  // skip that number + 3
		  for (int i = 0; i < linesToSkip+3; i++)
		  {
			  getline (tempFile,line);
		  }


		  // // read number of parameter values
		  getline (tempFile,line);
		  ss<<line;
		  linesToSkip;
		  ss>>linesToSkip;
		  ss=stringstream();

		  // skip that number + 2
		  for (int i = 0; i < linesToSkip+2; i++)
		  {
			  getline (tempFile,line);
		  }

		  // // read number of geometry entity indices
		  getline (tempFile,line);
		  ss<<line;
		  linesToSkip;
		  ss>>linesToSkip;
		  ss=stringstream();

		  // skip that number + 2
		  for (int i = 0; i < linesToSkip+2; i++)
		  {
			  getline (tempFile,line);
		  }

		  // // read number of up down pairs
		  getline (tempFile,line);
		  ss<<line;
		  linesToSkip;
		  ss>>linesToSkip;
		  ss=stringstream();

		  // skip that number + 8
		  for (int i = 0; i < linesToSkip+7; i++)
		  {

			  getline (tempFile,line);
				
		  }

		  //read  type of element 
		  getline (tempFile,line);
		  ss<<line;
		   ss>>nodesPerElement;
		  ss=stringstream();


		  // // read number of elements


		  getline (tempFile,line);
		  ss<<line;
		  int numElements;
		  ss>>numElements;
		  ss=stringstream();

		  meshConnectivity.resize(numElements*nodesPerElement); // each tri element has 6 nodes

		  // read mesh connectivity
		  getline (tempFile,line);
		  for (int i = 0; i < numElements; i++)
		  {
			   getline (tempFile,line);
			   ss<<line;
			   for (int j = 0; j < nodesPerElement; j++)
			   {
				   ss>>meshConnectivity[i*nodesPerElement+j ]; // again there are 6 nodes per element

			   }
			    
			   
		  	  ss=stringstream();
		  }

       tempFile.close();
	   cout<<" num nodes " << numNodes<<endl;
	   cout<< " num corners " <<numCorners <<endl;
	   cout<<" num elements " << numElements <<endl;
	  }
	else{

	cout<< " the file with mesh info could not be opened" <<endl;
	}
}
