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
#define TIME_STEP_DEBUG 1000
#include <math.h>
#include <float.h>
using namespace std;
using namespace Eigen;
//#define PRINT_DEBUG_INFO
//#define PRINT_DEBUG_INFO_GLOBAL_TIMING
//#define PRINT_DEBUG_INFO_NEIGHBORS_TIMING

const double pi =M_PI;
// imitating daniel's program structure Hinge

struct Hinge
{
public:
	int fiber_ID;
	int numberOfSegments;
	Vector3d position; 
	double length;
	bool isStationary;
	Vector3d velocity;
	Vector3d vorticity;
	Vector3d fluid_velocity_sum;
	Vector3d torque;
	double angle;
	int numberOfBeads;
	//Vector3d excludedVolumeForce;
	Vector3d sumOmegaFluid;
	Vector3d r_sum;
	Matrix3d r_product_sum;
	Matrix3d r_times_u_sum;
	Vector3d r;
	double averageviscosity;
	Vector3d exluded_volume_force;
	Vector3d excluded_volume_torque;
	int elementID;

};

