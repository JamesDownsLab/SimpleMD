// SimpleMD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Engine.h"
#include "setup.h"
#define FMT_HEADER_ONLY
#include <fmt/format.h>





int main()
{
    // setup
    //DimpleTest(1e-3, 10, 1e-6, 1e-9);
    double dimple_rad = 0.2e-3;
    double D = 1e-8;
    double Dend = 1e-8;
    HexGrid(4.8e-3, 5, 5, 1e-6, D, 5, 5, 1.0, dimple_rad, 100);
    fs::path savepath = "C:\\Users\\james\\Data\\output.dump";
    std::cout << savepath.string() << "\n";
    ProgramOptions options{
        savepath,
        Integrator::VerletPosition,
        Optimiser::Lattice,
        /* save_interval */ 1000,
        /*save_on_collision*/ false,
        /*random_force*/ true,
        /*seed*/ 2
    };

    Engine engine("initial.random", options);

    int smallStep = 1e4;
    int NsmallSteps = 10;
    double dD = (D - Dend) / (NsmallSteps - 1);

    int i{ 0 };
    while (i < NsmallSteps*smallStep) {
        for (int j{ 0 }; j < smallStep; j++) engine.step();
        //D -= dD;
        //engine.set_noise(D);
        i += smallStep;
        std::cout << "Step: " << i << " collisions: " << engine.collisions() << " Energy: " << engine.total_kinetic_energy() << " Force: " << engine.total_force() << "\n";
    }
}

//int main()
//{
//	double timestep = 1e-6;
//	double dimple_rad = 1e-3;
//	std::vector<double> dimple_ks{ 1, 10, 100 };
//	std::vector<double> Ds{ 1e-9, 1e-8, 1e-7 };
//	fs::path savepath = "C:\\Users\\james\\Data\\output.dump";
//	    ProgramOptions options{
//	        savepath,
//	        Integrator::VerletPosition,
//	        Optimiser::Lattice,
//	        /* save_interval */ 1000,
//	        /*save_on_collision*/ false,
//	        /*random_force*/ true,
//	        /*seed*/ 2
//	    };
//	for (double k : dimple_ks) {
//		for (double D : Ds) {
//			DimpleTest(dimple_rad, k, timestep, D);
//			std::string fname = fmt::format("Dimple k({0}) D({1}).dump", k, D);
//			savepath.replace_filename(fname);
//			std::cout << savepath << std::endl;
//			options.savepath = savepath;
//			Engine engine("initial.random", options);
//			int i{ 0 };
//			while (i < 1e6) {
//				for (int j{ 0 }; j < 1e4; j++) engine.step();
//				i += 1e4;
//				std::cout << "Step: " << i << "\r";
//			}
//			std::cout << "\n";
//		}
//	}
//}