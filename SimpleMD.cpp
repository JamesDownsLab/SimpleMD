// SimpleMD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Engine.h"
#include "setup.h"



int main()
{
	// setup
	Balls100();

	Engine engine("initial.random", 1e10, Integrator::VerletPosition, true);

    int i{ 0 };
    while (engine.collisions() < 20000) {
        for (int j{ 0 }; j < 1000; j++) engine.step();
        i += 1000;
        std::cout << "Step: " << i << " collisions: " << engine.collisions() << " Energy: " << engine.total_kinetic_energy() << "\r";
    }
}