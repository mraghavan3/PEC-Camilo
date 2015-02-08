#include "fiber.h"


//When a fiber breaks, a new fiber is created and the hinges of the old fiber are assigned to the 
// new fiber. 
fiber fiber::splitFiber(int hingeID){

	fiber newFiber; 

	if (hingeID<this->numberOfHinges)
	{
		for (int i = hingeID; i < this->numberOfHinges; i++)
		{
			newFiber.hinges.push_back(this->hinges[i]); // Copies the hinges that will be removed to the new created fiber

		}

		newFiber.numberOfHinges= this-> numberOfHinges - hingeID;
		newFiber.isBroken =true;
		this->hinges.resize(hingeID+1); // destroys the hinges that were created
		this->numberOfHinges = hingeID+1;
		this->isBroken =true;
	}
	else
	{
		cout<< " Problem splitting fiber";

	}

	
	return newFiber;
}
