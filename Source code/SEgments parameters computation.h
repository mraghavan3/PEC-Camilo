#pragma once
#include "classes and structures.h"
#include "fiber.h"


void fibers_parameters_computation(vector<fiber> &fibers);
void fibers_parameters_computation(vector<fiber> & fibers, double beadRadius);
void fibers_parameters_computation_mesh(vector<fiber> & fibers, double beadRadius, vector<double > &meshPoints, vector<int> &meshConnectivity, vector<double> &meshVelocities, vector<double> &meshVorticity, int numElements, int numNodes, int nodesPerElement );
