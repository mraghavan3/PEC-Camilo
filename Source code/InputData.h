#pragma once
#include <map>
#include <vector>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string> 
#include <sstream> 
#include <Eigen\Dense>
#define _USE_MATH_DEFINES
#include <math.h>
#include "classes and structures.h"
#include "fiber.h"

using namespace std;
using namespace Eigen;


int readBorder(std::vector<double>  & border_coords, string filename);
int readVelocities(std::vector<double>  & flowField, string filename);
int readViscosity(std::vector<double>  & flowField, string filename);
int readVorticity(std::vector<double>  & flowField, string filename);
int readFlowFieldData(std::vector<double>  & flowField, string filename);
int readInitialFibersPositions(vector<fiber> &fibers, vector<Hinge> &hinges);
int printPositions(vector<fiber> fibers, ofstream &outputFile, int frames);
void readVelFile2D(string Filename, vector<double> &velValues, vector<double> &vorticity  ); // reads data from comsol file
void readMeshFileTri(string Filename, vector<int> &meshConnectivity, vector<double> &meshPoints, vector<int> & meshCorners , int & nodesPerElement  );// reads mesh data from comsol file
void generateCornerFibers(vector<int > corners,vector<fiber>&fibers,vector<double> meshPoints, double depth, int numNodes  );

