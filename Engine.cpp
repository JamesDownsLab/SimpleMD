#include "Engine.h"

void Engine::step()
{
	collision = false;
	std::for_each(particles.begin(), particles.end(), [&](auto& p) {p.reset_contact(); });
	
	integrate();
	auto result = std::find_if(particles.begin(), particles.end(), [](auto& p) { return p.contact(); });
	collision = result == std::end(particles) ? false : true;
	std::for_each(particles.begin(), particles.end(), [](auto& p) {p.update_collisions(); });
	check_dump();
}

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

void Engine::init_system(const char* fname)
{
	std::ifstream fparticle{ fname };
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
		else if (type == "#nstep:") {
			fparticle >> nstep;
			fparticle.ignore(100, '\n');
			std::cout << "nstep: " << nstep << std::endl;
		}
		else if (type == "#timestep:") {
			fparticle >> timestep;
			fparticle.ignore(100, '\n');
			std::cout << "timestep: " << timestep << std::endl;
		}
		else if (type == "#nprint:") {
			fparticle >> nprint;
			fparticle.ignore(100, '\n');
			std::cout << "nprint: " << nprint << std::endl;
		}
		else if (type == "#nenergy:") {
			fparticle >> nenergy;
			fparticle.ignore(100, '\n');
			std::cout << "nenergy: " << nenergy << std::endl;
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
		else {
			std::cerr << "init: unknown global property: " << type << std::endl;
			abort();
		}
	}
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

void Engine::make_forces()
{
	for (unsigned int i{ 0 }; i < no_of_particles; i++) {
		for (unsigned int k{ i+1 }; k < no_of_particles; k++) {
			force(particles[i], particles[k], lx, ly);
		}
	}
}

void Engine::integrate()
{
	// Set forces to zero and move particles that belong to moving boundaries
	std::for_each(particles.begin(), particles.end(),
		[&](Particle& p) {
			if (p.pstate() == ParticleState::Free) {
				p.set_force_to_zero();
				if (integrator == Integrator::Gear) p.gear_predict(timestep);
			}
			else {
				p.boundary_conditions(timestep, Time);
			}
		});

	// Calculate all the forces between particles
	make_forces();

	//if (apply_potential) {

	//	G = Vector(p0 * cos(dis(gen)), p0 * sin(dis(gen)));
	//}
	//else G = null;

	// Update the positions of all the particles 
	std::for_each(particles.begin(), particles.end(),
		[&](Particle& p) {
			if (p.pstate() == ParticleState::Free) {
				switch (integrator)
				{
				case Integrator::VerletPosition:
					G = Vector(dis(gen), dis(gen));
					p.position_verlet(timestep, G);
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

void Engine::check_dump()
{
	if (save != save_interval) {
		save++;
	}
	else {
		//if (!collision) {
		save = 1;
		dump();
		//}
	}
	_collisions = collisions();
	if (_collisions != _last_collisions) {
		dump();
	}
	_last_collisions = _collisions;
}

void Engine::dump()
{
	std::fprintf(f1, "ITEM: TIMESTEP\n%d\n", int(Time / timestep));
	std::fprintf(f1, "ITEM: TIME\n\%.8f\n", Time);
	std::fprintf(f1, "ITEM: BOX BOUNDS pp pp f\n%.4f %.4f\n%.4f %.4f\n-0.002 0.002\n", x_0, x_0 + lx, y_0, y_0 + ly);
	std::fprintf(f1, "ITEM: NUMBER OF ATOMS\n%d\n", no_of_particles);
	std::fprintf(f1, "ITEM: KINETIC ENERGY\n%.3e\n", total_kinetic_energy());
	std::fprintf(f1, "ITEM: COLLISION\n%d\n", collision);
	std::fprintf(f1, "ITEM: ATOMS x y z vx vy radius type contact\n");
	for (Particle& p : particles) {
		std::fprintf(f1, "%.9f %.9f %.5f %.5f %.5f %.5f %d %d\n", p.x(), p.y(), 0.0, p.vx(), p.vy(), p.r(), p.pstate(), p.contact());
	}
}
