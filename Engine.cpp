#include "Engine.h"

///////////////////////////////////////////////////////////////////////////////
/// Initialisation
///////////////////////////////////////////////////////////////////////////////

void Engine::init_system(const char* fname)
{
  std::ifstream fparticle{ fname };
  // Read the system properties
  while (fparticle.peek() == '#') {
    std::string type;
    fparticle >> type;
    if (type == "#gravity:") {
      fparticle >> G.x() >> G.y();
      fparticle.ignore(100, '\n');
      std::cout << "gravity: " << G << std::endl;
    }
    else if (type == "#Time:") {
      fparticle >> Time;
      fparticle.ignore(100, '\n');
      std::cout << "Time: " << Time << std::endl;
    }
    else if (type == "#timestep:") {
      fparticle >> timestep;
      fparticle.ignore(100, '\n');
      std::cout << "timestep: " << timestep << std::endl;
    }
    else if (type == "#lx:") {
      fparticle >> lx;
      fparticle.ignore(100, '\n');
      std::cout << "lx: " << lx << std::endl;
    }
    else if (type == "#ly:") {
      fparticle >> ly;
      fparticle.ignore(100, '\n');
      std::cout << "ly: " << ly << std::endl;
    }
    else if (type == "#x_0:") {
      fparticle >> x_0;
      fparticle.ignore(100, '\n');
      std::cout << "x_0: " << x_0 << std::endl;
    }
    else if (type == "#y_0:") {
      fparticle >> y_0;
      fparticle.ignore(100, '\n');
      std::cout << "y_0: " << y_0 << std::endl;
    }
    else if (type == "#D:") {
      fparticle >> noise_strength;
      fparticle.ignore(100, '\n');
      std::cout << "D: " << noise_strength << std::endl;
    }
    else if (type == "#r_d:") {
      fparticle >> dimple_rad;
      fparticle.ignore(100, '\n');
      std::cout << "dimple rad: " << dimple_rad << std::endl;
    }
    else if (type == "#k_d:") {
      fparticle >> dimple_k;
      fparticle.ignore(100, '\n');
      std::cout << "dimple k: " << dimple_k << std::endl;
    }
    else if (type == "#DL:") {
      fparticle >> dimple_spacing;
      fparticle.ignore(100, '\n');
      std::cout << "dimple spacing: " << dimple_spacing << std::endl;
    }
    else {
      std::cerr << "init: unknown global property: " << type << std::endl;
      abort();
    }
  }

  // Read the dimples
  while (fparticle.peek() == '%') {
    std::string type;
    fparticle >> type;
    double x, y;
    if (type == "%dimple:") {
      fparticle >> x >> y;
      fparticle.ignore(100, '\n');
      dimples.emplace_back(x, y);
    }
  }

  // Read the particles
  while (fparticle) {
    Particle pp;
    fparticle >> pp;
    if (fparticle) {
      particles.push_back(pp);
    }
  }
  no_of_particles = particles.size();
  std::cout << no_of_particles << " particles read" << std::endl;

  std::for_each(particles.begin(), particles.end(), [&](auto& p) {p.init_rtd0_old(timestep); });

  dump();
}

///////////////////////////////////////////////////////////////////////////////
/// Main Steps
///////////////////////////////////////////////////////////////////////////////


Engine::Engine(const char* fname, ProgramOptions options) : _options{options}
{
  f1 = fopen(_options.savepath.string().c_str(), "w");
  rng.seed(options.seed);
  init_system(fname);
  if (options.optimiser == Optimiser::LinkCell) { 
    init_link_cell_algorithm(); 
  }
  if (_options.optimiser == Optimiser::Lattice) { 
    init_lattice_algorithm(); 
  }
}

void Engine::step()
{
  // Check whether the optimiser needs updating
  if (_options.optimiser == Optimiser::LinkCell) { make_link_cell(); }
  else if (_options.optimiser == Optimiser::Lattice) {
    if (ilist_needs_update()) { make_ilist(); }
  }
  
  // Collision counting
  collision = false;
  std::for_each(particles.begin(), particles.end(), [&](auto& p) {p.reset_contact(); });
  
  integrate();

  // More collision counting
  auto result = std::find_if(particles.begin(), particles.end(), [](auto& p) { return p.contact(); });
  collision = result == std::end(particles) ? false : true;
  std::for_each(particles.begin(), particles.end(), [](auto& p) {p.update_collisions(); });


  check_dump();
}

void Engine::integrate()
{
  // Set forces to zero and move particles that belong to moving boundaries
  std::for_each(particles.begin(), particles.end(),
    [&](Particle& p) {
      if (p.pstate() == ParticleState::Free) {
        p.set_force_to_zero();
        if (_options.integrator == Integrator::Gear) p.gear_predict(timestep);
      }
      else {
        p.boundary_conditions(timestep, Time);
      }
    });

  // Calculate all the forces between particles
  make_forces();
  if (_options.random_force == true) {
    make_random_forces();
    correct_random_forces();
  }
  calculate_dimple_force();


  // Update the positions of all the particles 
  std::for_each(particles.begin(), particles.end(),
    [&](Particle& p) {
      if (p.pstate() == ParticleState::Free) {
        switch (_options.integrator)
        {
        case Integrator::VerletPosition:
          p.position_verlet(timestep, lx, ly, G);
          break;
        case Integrator::VerletVelocity:
          p.velocity_verlet(timestep, G);
          break;
        case Integrator::Gear:
          p.gear_correct(timestep, G);
          break;
        default:
          break;
        }
      }
    });

  // Apply periodic boundary conditions
  std::for_each(particles.begin(), particles.end(),
    [&](Particle& p) {p.periodic_bc(x_0, y_0, lx, ly); });
  Time += timestep;
}

///////////////////////////////////////////////////////////////////////////////
/// Calculating Forces
///////////////////////////////////////////////////////////////////////////////

void Engine::make_forces()
{
  if (_options.optimiser == Optimiser::LinkCell) {
    // for all cells
    for (unsigned int ix{ 0 }; ix < linkCell.size(); ix++) {
      for (unsigned int iy{ 0 }; iy < linkCell[ix].size(); iy++) {
        for (unsigned int j{ 0 }; j < linkCell[ix][iy].size(); j++) {
          int pj = linkCell[ix][iy][j];

          // force between all particles in the same cell
          for (unsigned int k{ j + 1 }; k < linkCell[ix][iy].size(); k++) {
            int pk = linkCell[ix][iy][k];
            force(particles[pj], particles[pk], lx, ly);
          }

          // force between particles in neighbouring cells
          for (unsigned int n{ 0 }; n < neighbours[ix][iy].size(); n++) {
            int iix = neighbours[ix][iy][n].first;
            int iiy = neighbours[ix][iy][n].second;

            for (unsigned int k{ 0 }; k < linkCell[iix][iiy].size(); k++) {
              int pk = linkCell[iix][iiy][k];
              force(particles[pj], particles[pk], lx, ly);
            }
          }
        }
      }
    }
  }
  else if (_options.optimiser == Optimiser::Lattice) {
    // Loop over the partners list for each particle
    for (unsigned int i{ 0 }; i < no_of_particles; i++) {
      for (unsigned int k{ 0 }; k < partners[i].size(); k++) {
        int pk = partners[i][k];
        force(particles[i], particles[pk], lx, ly);
      }
    }
  }
  else {
    // Loop over all pairs of particles
    for (unsigned int i{ 0 }; i < no_of_particles; i++) {
      for (unsigned int k{ i + 1 }; k < no_of_particles; k++) {
        force(particles[i], particles[k], lx, ly);
      }
    }
  }

}

void Engine::make_random_forces()
{
  // Add a random force to each particle using David Bray's method
  for (auto& p : particles) {
    double a1 = a_dis();
    double a2 = a_dis();
    double K = sqrt(-4 * log(a1) * noise_strength / timestep);
    Vector f{ K * cos(2 * PI * a2), K * sin(2 * PI * a2) };
    p.set_random_force(f);
  }
}

void Engine::correct_random_forces()
{
  // Calculate the net random force
  Vector total_force{ null };
  for (auto& p : particles) {
    total_force += p.random_force();
  }

  // Subtract a share of this force from each particle
  total_force *= (-1.0 / no_of_particles);
  for (auto& p : particles) {
    p.correct_random_force(total_force);
  }
}

void Engine::calculate_dimple_force()
{
  double dlx = sqrt(3) * dimple_spacing;
  double dly = dimple_spacing;
  for (auto& p : particles)
  {
    double dx = fmod(p.x(), dlx);
    double dy = fmod(p.y(), dly);
    if (dx < dlx / 2) {
      if (dy < dly / 2) {
        double r2_a = dx * dx + dy * dy;
        double r2_b = (dlx / 2 - dx) * (dlx / 2 - dx) + (dly / 2 - dy) * (dly / 2 - dy);
        if (r2_a < r2_b) {
          dx = -dx;
          dy = -dy;
        }
        else {
          dx = (dlx / 2 - dx);
          dy = (dly / 2 - dy);
        }
      }
      else {
        double r2_a = dx * dx + (dly - dy) * (dly - dy);
        double r2_b = (dlx - dx) * (dlx - dx) + (dy - dly / 2) * (dy - dly / 2);
        if (r2_a < r2_b) {
          dx = -dx;
          dy = (dly - dy);
        }
        else {
          dx = (dlx / 2 - dx);
          dy = -(dy - dly / 2);
        }
      }
    }
    else {
      if (dy < dly / 2) {
        double r2_a = (dlx - dx) * (dlx - dx) + dy * dy;
        double r2_b = (dlx / 2 - dx) * (dlx / 2 - dx) + (dly / 2 - dy) * (dly / 2 - dy);
        if (r2_a < r2_b) {
          dx = dlx - dx;
          dy = -dy;
        }
        else {
          dx = -(dx - dlx / 2);
          dy = dly / 2 - dy;
        }
      }
      else {
        double r2_a = (dlx - dx) * (dlx - dx) + (dly - dy) * (dly - dy);
        double r2_b = (dlx / 2 - dx) * (dlx / 2 - dx) + (dy - dly / 2) * (dy - dly / 2);
        if (r2_a < r2_b) {
          dx = dlx - dx;
          dy = dly - dy;
        }
        else {
          dx = -(dx - dlx / 2);
          dy = -(dy - dly / 2);
        }
      }
    }
    if (abs(dx) < p.x() + dimple_rad && abs(dy) < p.y() + dimple_rad) {
      double rr = sqrt(dx * dx + dy * dy);
      if (rr > 0) {
        double xi = p.r() + dimple_rad - rr;
        double rr_rez = 1 / rr;
        double ex = dx * rr_rez;
        double ey = dy * rr_rez;
        // Force
        double elastic_force = dimple_k * xi;
        double fn = elastic_force;
        //if (fn > 0) fn = 0;
        if (p.pstate() == ParticleState::Free) {
          p.add_dimple_force(Vector(fn * ex, fn * ey));
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/// Lattice Method 
///////////////////////////////////////////////////////////////////////////////

void Engine::init_lattice_algorithm()
{
  // Calcualte gk, the size of the lattice sites
  rmin = particles[0].r();
  rmax = particles[0].r();
  for (auto& p : particles) {
    if (p.r() < rmin) rmin = p.r();
    if (p.r() > rmax) rmax = p.r();
  }
  gk = sqrt(2) * rmin;

  // Calculate gm, the number of lattice sites that the largest particle covers
  gm = int(2 * rmax / gk) + 1;

  // Calculate the number of lattice sites in each dimension
  Nx = int(lx / gk) + 1;
  Ny = int(ly / gk + 1);
  partners.resize(no_of_particles);
  pindex.resize(Nx);
  for (auto& p : pindex) {
    p.resize(Ny);
  }
  clear_pindex();
  make_ilist();
}

void Engine::make_ilist()
{
  // For each particle, add it to the pindex lattice site
  // and clear its partners
  for (unsigned int i{ 0 }; i < no_of_particles; i++) {
    double x = particles[i].x();
    double y = particles[i].y();
    int ix = int((x - x_0) / gk);
    int iy = int((y - y_0) / gk);
    pindex[ix][iy] = i;
    partners[i].clear();
  }

  // Generate the partners list for each particle
  for (unsigned int i{ 0 }; i < no_of_particles; i++) {
    double x = particles[i].x();
    double y = particles[i].y();
    if ((x >= x_0) && (x < x_0 + lx) && (y >= y_0) && (y < y_0 + ly)) {
      int ix = int((x - x_0) / gk);
      int iy = int((y - y_0) / gk);
      // Check the adjacent gm lattice sites for particles
      for (int dx = -gm; dx <= gm; dx++) {
        for (int dy = -gm; dy <= gm; dy++) {
          int iix = (ix + dx + Nx) % Nx;
          int iiy = (iy + dy + Ny) % Ny;
          int k = pindex[iix][iiy];
          // Only record the particle once
          if (k > (int)i) {
            partners[i].push_back(k);
          }
        }
      }
    }
  }
}

bool Engine::ilist_needs_update() {
  // If the pindex is the same as last time then it doesn't need updating
  for (unsigned int i{ 0 }; i < no_of_particles; i++) {
    double x = particles[i].x();
    double y = particles[i].y();
    int ix = int((x - x_0) / gk);
    int iy = int((y - y_0) / gk);
    if (pindex[ix][iy] != i) {
      clear_pindex();
      return true;
    }
  }
  return false;
}

void Engine::clear_pindex()
{
  for (auto& p : pindex) {
    for (auto& q : p) {
      q = -1;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/// Link Cell
///////////////////////////////////////////////////////////////////////////////

void Engine::make_link_cell()
{
  // Lists are emptied
  for (unsigned int ix{ 0 }; ix < linkCell.size(); ix++) {
    for (unsigned int iy{ 0 }; iy < linkCell[ix].size(); iy++) {
      linkCell[ix][iy].clear();
    }
  }

  // Assign each particle a box which contains its center.
  for (int i{ 0 }; i < no_of_particles; i++) {
    int ix = int(nx * (particles[i].x() - x_0) / lx);
    int iy = int(ny * (particles[i].y() - y_0) / ly);
    if ((ix >= 0) && (ix < nx) && (iy >= 0) && (iy < ny)) {
      linkCell[ix][iy].push_back(i);
    }
    else {
      std::cerr << "Particle " << i << " outside simulation area\n";
      exit(0);
    }
  }
}

bool Engine::is_valid_neighbour(int ix, int iy, int iix, int iiy)
{
  // check whether boxes have to be recorded.
  if ((iix == (ix - 1 + nx) % nx) && (iiy == (iy + 1 + ny) % ny)) return true;
  if ((iix == (ix + nx) % nx) && (iiy == ((iy + 1 + ny) % ny))) return true;
  if ((iix == (ix + 1 + nx) % nx) && (iiy == (iy + 1 + ny) % ny)) return true;
  if ((iix == (ix + 1 + nx) % nx) && (iiy == (iy + ny) % ny)) return true;
  return false;
}

void Engine::init_neighbours()
{
  int iix, iiy;

  neighbours.resize(nx);
  for (int ix{ 0 }; ix < nx; ix++) {
    neighbours[ix].resize(ny);
  }

  for (int ix{ 0 }; ix < nx; ix++) {
    for (int iy{ 0 }; iy < ny; iy++) {
      for (int dx{ -1 }; dx <= 1; dx++) {
        for (int dy{ -1 }; dy < 1; dy++) {
          iix = (ix + dx + nx) % nx;
          iiy = (iy + dy + ny) % ny;
          if (is_valid_neighbour(ix, iy, iix, iiy)) {
            neighbours[ix][iy].push_back(std::pair<int, int>(iix, iiy));
          }
        }
      }
    }
  }
}

void Engine::init_link_cell_algorithm()
{
  linkCell.resize(nx);
  for (unsigned int ix{ 0 }; ix < linkCell.size(); ix++) {
    linkCell[ix].resize(ny);
  }
  init_neighbours();
  make_link_cell();
}

///////////////////////////////////////////////////////////////////////////////
/// Dumping
///////////////////////////////////////////////////////////////////////////////

void Engine::check_dump()
{
  if (save != _options.save_interval) {
    save++;
  }
  else {
    //if (!collision) {
    save = 1;
    dump();
    //}
  }
  if (_options.save_on_collision == true) {
    _collisions = collisions();
    if (_collisions != _last_collisions) {
      dump();
    }
    _last_collisions = _collisions;
  }
}

void Engine::dump()
{
  std::fprintf(f1, "ITEM: TIMESTEP\n%d\n", int(Time / timestep));
  std::fprintf(f1, "ITEM: TIME\n%.8f\n", Time);
  std::fprintf(f1, "ITEM: BOX BOUNDS pp pp f\n%.4f %.4f\n%.4f %.4f\n-0.002 0.002\n", x_0, x_0 + lx, y_0, y_0 + ly);
  std::fprintf(f1, "ITEM: NUMBER OF ATOMS\n%d\n", no_of_particles + dimples.size());
  std::fprintf(f1, "ITEM: KINETIC ENERGY\n%.3e\n", total_kinetic_energy());
  std::fprintf(f1, "ITEM: COLLISION\n%d\n", collision);
  std::fprintf(f1, "ITEM: ATOMS x y z vx vy radius type contact\n");
  for (Particle& p : particles) {
    std::fprintf(f1, "%.9f %.9f %.5f %.5f %.5f %.5f %d %d\n", p.x(), p.y(), 0.0, p.vx(), p.vy(), p.r(), p.pstate(), p.contact());
  }
  for (Vector& d : dimples) {
    std::fprintf(f1, "%.9f %.9f %.5f %.5f %.5f %.5f %d %d\n", d.x(), d.y(), 0.0, 0.0, 0.0, dimple_rad, 2, 0);
  }
}

///////////////////////////////////////////////////////////////////////////////
/// Extra Calculations
///////////////////////////////////////////////////////////////////////////////

int Engine::collisions()
{
  int sum{ 0 };
  std::for_each(particles.begin(), particles.end(), [&](auto& p) {sum += p.collisions(); });
  return sum;
}

double Engine::total_kinetic_energy()
{
  double energy{ 0 };
  for (auto& p : particles) {
    energy += p.kinetic_energy();
  }
  return energy;
}

double Engine::total_force()
{
  double f{ 0 };
  for (auto& p : particles) {
    f += p.force().length();
  }
  return f;
}
