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


enum coordinates{ X, Y, Z};
// " global variables " these should be handled better

//  Stores the information of the border of the domain
//  x1 y1 z1 x2 y2 z2 for each segment width= 6
	std::vector<double>  border_coords;

// Stores information about the flow field
//  x y z Vx Vy Vz vorticity viscosity width=8
	std::vector<double> flowField;

// input variables

	double viscosity;
	double fiber_radius;
	double E_young;
	double dt;
	double L_segment;
	double max_alpha;
	double gamma_dot;
	double critical_angle;
	int number_fibers;
	int number_integrations;
	int writing_frequency;
	int check_breakage_frequency;

	enum variablesNamesValues  // this is for reading the data in 
	{
		unkownVariable,fiber_radius_name, viscosity_name , 
		E_young_name, dt_name, number_integrations_name,
		gamma_dot_name, writing_frequency_name, L_seg_name,
		critical_angle_name, check_breakage_frequency_name,
		
	};
	// Map to associate the strings with the enum values avobe
static std::map<std::string, variablesNamesValues> s_mapVariablesValues;





vector<Hinge> hinges;