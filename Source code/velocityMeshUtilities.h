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




int find_tri_neighbors(int numTris, int nodesPerTri, vector<int> &tri_points, vector<int> & triNeighbors);
void find_Initial_tetID(vector<double>&velMeshPoints,vector<int>&velMeshConnectivity,vector<fiber> &fibers, int numElements, int numnodes,int meshType,int nodesPerElement );
void findBoundaryEdges(int numElements, int numFacesPerElement, vector<int>& elementNeighbors, vector<int > & boundaryFaces);
void checkHingesElements (const vector<double>&velMeshPoints,const vector<int>&velMeshConnectivity, const vector<int> meshNeighbors,vector<fiber> &fibers, int numElements, int numnodes,int meshType,int nodesPerElement );