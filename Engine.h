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

class Engine
{
public:
	Engine(const char* fname, int save_int, Integrator inter, bool potential = false) : save_interval{ save_int }, integrator{ inter }, apply_potential{ potential } {
		init_system(fname);
	};

	void step();
	int collisions();
	double total_kinetic_energy();

private:

	void init_system(const char* fname);
	void make_forces();
	void integrate();
	void check_dump();
	void dump();



	double Time{ 0 };
	double timestep;
	bool collision{ false };
	int _last_collisions{ 0 };
	int _collisions{ 0 };
	bool apply_potential{ false };


	int nstep;
	int nprint;
	int nenergy;
	int save{ 1 };
	int save_interval;

	double lx;
	double ly;
	double x_0;
	double y_0;


	const double PI = 3.1415926;
	std::random_device rd;
	std::mt19937 gen{ rd() };
	double p0{ 100.0 };
	std::normal_distribution<double> dis{ 0.0, p0 };

	Integrator integrator;

	Vector G;

	unsigned int no_of_particles;

	FILE* f1 = std::fopen("C:/Users/james/Data/output.dump", "w");

	std::vector<Particle> particles;
};

