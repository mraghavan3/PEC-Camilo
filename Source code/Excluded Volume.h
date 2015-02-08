#pragma once
#include "classes and structures.h"
#include "Broadphase.h"
#include "fiber.h"
#include "utilities.h"

void refineNeighbors(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  possibleCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices);
void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius);
void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  possibleCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices  ); 
void calculateFiberInteractions(vector<fiber> & fibers, double beadRadius, thrust::host_vector<long long>  potentialCollisions, thrust::host_vector<int>  fiberIndices, thrust::host_vector<int>  segmentIndices , int timeStep );
void interactionsWallsTri(vector<fiber> & fibers , int numNodes,vector<double> & meshPoints , vector<int> &meshConnectivity , vector<int> meshNeighbors , double fiberRadius , double excludedVolumeForce, double threshold, int nodesPerElements);
