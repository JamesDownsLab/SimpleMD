// SimpleMD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Engine.h"
#include "setup.h"





int main()
{
    // setup
    DimpleTest();

    ProgramOptions options{
        Integrator::VerletPosition,
        Optimiser::Lattice,
        /* save_interval */ 1000,
        /*save_on_collision*/ false,
        /*random_force*/ true,
        /*seed*/ 2
    };

    Engine engine("initial.random", options);

    int i{ 0 };
    while (i < 10000000) {
        for (int j{ 0 }; j < 1e3; j++) engine.step();
        i += 1000;
        std::cout << "Step: " << i << " collisions: " << engine.collisions() << " Energy: " << engine.total_kinetic_energy() << " Force: " << engine.total_force() << "\n";
    }
}
