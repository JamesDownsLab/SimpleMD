#pragma once

#include <iostream>
#include "Vector.h"
#include "ParticleState.h"

const double PI = 3.1415926;

inline double normalize(double dx, double L) {
	while (dx < -L / 2) dx += L;
	while (dx >= L / 2) dx -= L;
	return dx;
}
class Particle
{
	friend double Distance(const Particle& p1, const Particle& p2, double lx, double ly);
	friend void force(Particle& p1, Particle& p2, double lx, double ly);
	friend std::istream& operator >> (std::istream& is, Particle& p);
	friend std::ostream& operator << (std::ostream& os, const Particle& p);

public:
	// Initialisers
	Particle() : rtd0(null), rtd1(null), rtd2(null), rtd3(null), rtd4(null) {}
	Particle(const Particle& rhs) = default;
	Particle(Particle && rhs) = default;
	Particle& operator=(const Particle & rhs) = default;
	Particle& operator=(Particle && rhs) = default;
	virtual ~Particle() = default;

	// Getters
	double& r() { return _r; }
	double r() const { return _r; }
	double m() const { return _m; }
	ParticleState pstate() const { return _pstate; }
	const Vector& velocity() const { return rtd1; }
	double& x() { return rtd0.x(); }
	double x() const { return rtd0.x(); }
	double& y() { return rtd0.y(); }
	double y() const { return rtd0.y(); }
	double& vx() { return rtd1.x(); }
	double vx() const { return rtd1.x(); }
	double& vy() { return rtd1.y(); }
	double vy() const { return rtd1.y(); }
	bool& contact() { return _contact; }
	bool contact() const { return _contact; }
	double kinetic_energy() { return 0.5 * _m * (rtd1.x() * rtd1.x() + rtd1.y() * rtd1.y());}
	Vector& force() { return _force; }
	Vector force() const { return _force; }


	// Setters
	void reset_contact() { _contact = false; }
	void set_force_to_zero() { _force = null; _random_force = null; }
	void init_rtd0_old(double dt) { rtd0_old = rtd0 - rtd1 * dt; }

	Vector& random_force() { return _random_force; }
	Vector random_force() const { return _random_force; }
	void correct_random_force(Vector& correction) { _random_force += correction; }
	void add_force(const Vector& f) { _force += f; }
	void add_dimple_force(const Vector& f) { _dimple_force = f; }
	void set_random_force(Vector& f) { _random_force = f; }
	void periodic_bc(double x_0, double y_0, double lx, double ly);
	void boundary_conditions(double timestep, double Time);


	void gear_predict(double dt);
	void gear_correct(double dt, Vector G=null);
	void position_verlet(double dt, double lx, double ly, Vector G=null);
	void velocity_verlet(double dt, Vector G=null);

	void update_collisions();

	int& collisions() { return _collisions; }
	int collisions() const { return _collisions; }




private:
	Vector rtd0, rtd1, rtd2, rtd3, rtd4;
	Vector rtd0_old, rtd0_new;
	Vector _force;
	Vector _random_force{ null };
	Vector _dimple_force{ null };
	ParticleState _pstate{ ParticleState::Free };
	double _r, _m, _k, _epsilon, _gamma;
	bool _contact{ false }, _previous_contact{ false };
	int _collisions{ 0 };
};

