#pragma once

#include <algorithm>
#include <fstream>
#include <cstdio>
#include <vector>
#include <random>
#include <boost/random.hpp>
#include <math.h>
#include "Particle.h"
#include <filesystem>

const double SQRT3 = sqrt(3);

namespace fs = std::filesystem;
/**
 * The integration method used to move the particles.
 */
enum class Integrator {
  VerletPosition,
  VerletVelocity,
  Gear
};

/**
 * The optimsation methods for calculating the force.
 */
enum class Optimiser {
  /// Don't Optimise
  None,
  /// Link cell
  LinkCell,
  /// Lattice
  Lattice
};

struct ProgramOptions {
  fs::path savepath;
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
  /**
   * Initialises the Engine class.
   * 
   * \param fname File containing initialisation data
   * \param options struct containing various options for the program
   */
  Engine(const char* fname, ProgramOptions options);

  /// Iterates the simulation by one timestep.
  void step();

  /// Returns the total number of collisions at the current timestep.
  int collisions();

  /// Returns the total kinetic energy at the current timestep.
  double total_kinetic_energy();

  /// Returns the sum of the magnitudes of all forces from collisions.
  double total_force();

  /// Sets the constant D in the random force calculation.
  void set_noise(double s) {noise_strength = s;}


private:
  void init_system(const char* fname);

  /// Calculates the collisional forces between all particles.
  void make_forces();

  /// Generates the random force for each particle.
  void make_random_forces();

  /// Removes the net random force from all the particles to prevent drift.
  void correct_random_forces();

  /// Calculates forces and updates positions/velocities.
  void integrate();

  void check_dump();
  void dump();


  unsigned int no_of_particles;
  std::vector<Particle> particles;
  double Time{ 0 };
  double timestep;


  /////////////////////////////////////////////////////////////////////////////
  /// Collision Stuff
  /////////////////////////////////////////////////////////////////////////////

  /// True if the current timestep contains a collision
  bool collision{ false };
  /// The number of collisions on the previous timestep
  int _last_collisions{ 0 };
  /// The current number of collisions.
  int _collisions{ 0 };

  /////////////////////////////////////////////////////////////////////////////
  /// File Writing
  /////////////////////////////////////////////////////////////////////////////
  
  std::FILE* f1;
  int save{ 1 };

  /////////////////////////////////////////////////////////////////////////////
  /// System Properties
  /////////////////////////////////////////////////////////////////////////////=

  double lx;
  double ly;
  double x_0;
  double y_0;

  ProgramOptions _options;
  Vector G;

  /////////////////////////////////////////////////////////////////////////////
  /// Dimple Stuff
  /////////////////////////////////////////////////////////////////////////////
  void calculate_dimple_force();
  std::vector<Vector> dimples;
  double dimple_rad;
  double dimple_k, dimple_spacing;


  /////////////////////////////////////////////////////////////////////////////
  /// Random Force Stuff
  /////////////////////////////////////////////////////////////////////////////
  const double PI = 3.1415926;
  double noise_strength; // N^2s
  boost::mt19937 rng;
  boost::uniform_real<double> gen{ 0.0, 1.0 };
  boost::variate_generator <boost::mt19937&, boost::uniform_real<double>> a_dis{ rng, gen };

  /////////////////////////////////////////////////////////////////////////////
  /// Link Cell
  /////////////////////////////////////////////////////////////////////////////
  void make_link_cell();
  bool is_valid_neighbour(int ix, int iy, int iix, int iiy);
  void init_neighbours();
  void init_link_cell_algorithm();
  const int nx{ 5 }, ny{ 5 };
  std::vector<std::vector<std::vector<int>>> linkCell;
  std::vector<std::vector<std::vector<std::pair<int, int>>>> neighbours;
  

  /////////////////////////////////////////////////////////////////////////////
  /// Lattice
  /////////////////////////////////////////////////////////////////////////////
  
  /// 2D vector representing the lattice cells
  /// Value contains -1 if empty or the index of the particle.
  std::vector<std::vector<int>> pindex;

  /// Vector containing a list of neighbours for each particle.
  std::vector<std::vector<int>> partners;
  
  /// Updates pindex and partners.
  void make_ilist();

  /// Checks if pindex will change.
  bool ilist_needs_update();
  void clear_pindex();
  void init_lattice_algorithm();
  double rmin, rmax, gk;
  int gm, Nx, Ny;
  



};

