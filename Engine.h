#pragma once

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <random>
#include <math.h>

#include "Particle.h"

enum class Integrator {
	VerletPosition,
	VerletVelocity,
	Gear
};

enum class Optimiser {
	None,
	Verlet,
	LinkCell,
	LinkedList
};

struct ProgramOptions {
	Integrator integrator;
	Optimiser optimiser;
	double save_interval;
	bool save_on_collision;
	bool random_force;
	int seed;
};

class Engine
{
public:
	Engine(const char* fname, ProgramOptions options) : _options{ options }, gen{ options.seed } {
		init_system(fname);
		if (options.optimiser == Optimiser::LinkCell) { init_link_cell_algorithm(); }
		if (_options.optimiser == Optimiser::LinkedList) { init_lattice_algorithm(); }
	};

	void step();
	int collisions();
	double total_kinetic_energy();
	double total_force();


private:

	// Functions
	void init_system(const char* fname);
	void make_forces();
	void make_random_forces();
	void correct_random_forces();
	void integrate();
	void check_dump();
	void dump();


	// Data
	unsigned int no_of_particles;
	std::vector<Particle> particles;
	double Time{ 0 };
	double timestep;

	// Collision Recording
	bool collision{ false };
	int _last_collisions{ 0 };
	int _collisions{ 0 };

	// File Writing
	FILE* f1 = std::fopen("C:/Users/james/Data/output.dump", "w");
	int save{ 1 };

	// Temp File Writing
	FILE* f2 = std::fopen("C:/Users/james/Data/forces.txt", "w");
	int save2{ 1 };



	// Box Properties
	double lx;
	double ly;
	double x_0;
	double y_0;

	// Options
	ProgramOptions _options;
	Vector G;


	// Random Force Stuff
	const double PI = 3.1415926;
	std::mt19937 gen;
	double noise_strength; // N^2s
	std::uniform_real_distribution<double> a1_dis{ 0.0, 1.0 };
	std::uniform_real_distribution<double> a2_dis{ 0.0, 1.0 };

	// Link Cell Stuff
	void make_link_cell();
	bool is_valid_neighbour(int ix, int iy, int iix, int iiy);
	void init_neighbours();
	void init_link_cell_algorithm();
	const int nx{ 10 }, ny{ 10 };
	std::vector<std::vector<std::vector<int>>> linkCell;
	std::vector<std::vector<std::vector<std::pair<int, int>>>> neighbours;

	// Lattice Alg Stuff
	void make_ilist();
	void clear_pindex();
	void init_lattice_algorithm();
	double rmin, rmax, gk;
	int gm, Nx, Ny;
	std::vector<std::vector<int>> partners, pindex;



};

