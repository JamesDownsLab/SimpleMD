#pragma once

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <random>
#include <boost/random.hpp>
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
	Lattice
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
	Engine(const char* fname, ProgramOptions options) : 
		_options{ options }
	{
		rng.seed(options.seed);
		init_system(fname);
		if (options.optimiser == Optimiser::LinkCell) { init_link_cell_algorithm(); }
		if (_options.optimiser == Optimiser::Lattice) { init_lattice_algorithm(); }
		init_dimples();
	};

	void step();
	int collisions();
	double total_kinetic_energy();
	double total_force();
	void set_noise(double s) {noise_strength = s;}


private:

	// Functions
	bool init;
	void init_system(const char* fname);
	void init_dimples();
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



	// Box Properties
	double lx;
	double ly;
	double x_0;
	double y_0;

	// Options
	ProgramOptions _options;
	Vector G;

	// Dimple stuff
	void calculate_dimple_force();
	std::vector<Vector> dimples;
	double dimple_rad;
	std::vector<std::vector<std::vector<Vector>>> dimples_list;
	int nxd{ 10 }, nyd{ 10 };
	double dimple_k;


	// Random Force Stuff
	const double PI = 3.1415926;
	double noise_strength; // N^2s
	boost::mt19937 rng;
	boost::uniform_real<double> gen{ 0.0, 1.0 };
	boost::variate_generator <boost::mt19937&, boost::uniform_real<double>> a_dis{ rng, gen };

	// Link Cell Stuff
	void make_link_cell();
	bool is_valid_neighbour(int ix, int iy, int iix, int iiy);
	void init_neighbours();
	void init_link_cell_algorithm();
	const int nx{ 5 }, ny{ 5 };
	std::vector<std::vector<std::vector<int>>> linkCell;
	std::vector<std::vector<std::vector<std::pair<int, int>>>> neighbours;
	

	// Lattice Alg Stuff
	void make_ilist();
	bool ilist_needs_update();
	void clear_pindex();
	void init_lattice_algorithm();
	double rmin, rmax, gk;
	int gm, Nx, Ny;
	std::vector<std::vector<int>> partners, pindex;



};

