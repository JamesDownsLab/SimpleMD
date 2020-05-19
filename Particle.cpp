#include "Particle.h"

double Distance(const Particle& p1, const Particle& p2, double lx, double ly)
{
	double dx = normalize(p1.rtd0.x() - p2.rtd0.x(), lx);
	double dy = normalize(p1.rtd0.y() - p2.rtd0.y(), ly);
	return std::sqrt(dx * dx + dy * dy);
}

void force(Particle& p1, Particle& p2, double lx, double ly)
{
	double dx = normalize(p1.x() - p2.x(), lx);
	double dy = normalize(p1.y() - p2.y(), ly);
	if (abs(dx) < p1.r() + p2.r() && abs(dy) < p1.r() + p2.r()) {
		double rr = sqrt(dx * dx + dy * dy);
		double r1 = p1.r();
		double r2 = p2.r();

		// Overlap
		double xi = r1 + r2 - rr;

		if (xi > 0) { // If overlapping
			p1._contact = true;
			p2._contact = true;

			// Unit Vectors
			double rr_rez = 1 / rr;
			double ex = dx * rr_rez;
			double ey = dy * rr_rez;

			// Relative velocities
			double dvx = p1.vx() - p2.vx();
			double dvy = p1.vy() - p2.vy();

			// Overlap rate
			double xidot = -(ex * dvx + ey * dvy);

			// Average particle properties
			double k = 0.5 * (p1._k + p2._k);
			double gamma = 0.5 * (p1._gamma + p2._gamma);

			// Normal Force
			double elastic_force = k * xi;
			double dissipative_force = gamma * xidot;
			double fn = elastic_force + dissipative_force;

			// 
			if (fn < 0) fn = 0;

			if (p1.pstate() == ParticleState::Free) {
				p1.add_force(Vector(fn * ex, fn * ey));
			}
			if (p2.pstate() == ParticleState::Free) {
				p2.add_force(Vector(-fn * ex, -fn * ey));
			}
		}
	}
}

std::istream& operator>>(std::istream& is, Particle& p)
{
	is >> p.rtd0 >> p.rtd1
		>> p._r >> p._m >> p._pstate
		>> p._k >> p._epsilon
		>> p._force
		>> p.rtd2 >> p.rtd3 >> p.rtd4;
	p._gamma = abs(log(p._epsilon) / PI * sqrt((4 * (p._m*0.5) * p._k) / (1 + (log(p._epsilon) / PI) * (log(p._epsilon) / PI))));
	return is;
}

std::ostream& operator<<(std::ostream& os, const Particle& p)
{
	os << p.rtd0 << " " << p.rtd1 << " ";
	os << p._r << " " << p._m << " " << p._pstate << " ";
	os << p._k << " " << p._epsilon << " ";
	os << p._force << " ";
	os << p.rtd2 << " " << p.rtd3 << " " << p.rtd4 << "\n" << std::flush;
	return os;
}

void Particle::periodic_bc(double x_0, double y_0, double lx, double ly)
{
	while (rtd0.x() < x_0) rtd0.x() += lx;
	while (rtd0.x() > x_0 + lx) rtd0.x() -= lx;
	while (rtd0.y() < y_0) rtd0.y() += ly;
	while (rtd0.y() > y_0 + ly) rtd0.y() -= ly;
	//while (rtd0_old.x() < x_0) rtd0_old.x() += lx;
	//while (rtd0_old.x() > x_0 + lx) rtd0_old.x() -= lx;
	//while (rtd0_old.y() < y_0) rtd0_old.y() += ly;
	//while (rtd0_old.y() > y_0 + ly) rtd0_old.y() -= ly;
}

void Particle::boundary_conditions(double timestep, double Time)
{
	{
		switch (pstate())
		{
		case ParticleState::Free: break;
		case ParticleState::Confined: break;
		default:
			std::cerr << "pstate: " << pstate() << " not implemented\n";
			abort();
		}
	}
}

void Particle::gear_predict(double dt)
{
	double a1 = dt;
	double a2 = a1 * dt / 2;
	double a3 = a2 * dt / 3;
	double a4 = a3 * dt / 4;

	rtd0 += a1 * rtd1 + a2 * rtd2 + a3 * rtd3 + a4 * rtd4;
	rtd1 += a1 * rtd2 + a2 * rtd3 + a3 * rtd4;
	rtd2 += a1 * rtd3 + a2 * rtd4;
	rtd3 += a1 * rtd4;
}

void Particle::gear_correct(double dt, Vector G)
{
	static Vector accel, corr;
	double dtrez = 1 / dt;
	const double coeff0 = double(19) / double(90) * (dt * dt / double(2));
	const double coeff1 = double(3) / double(4) * (dt / double(2));
	const double coeff3 = double(1) / double(2) * (double(3) * dtrez);
	const double coeff4 = double(1) / double(12) * (double(12) * (dtrez * dtrez));

	accel = Vector((1 / _m) * (_force.x()+_random_force.x()) + G.x(),
		(1 / _m) * (_force.y()+_random_force.y()) + G.y());

	corr = accel - rtd2;
	rtd0 += coeff0 * corr;
	rtd1 += coeff1 * corr;
	rtd2 = accel;
	rtd3 += coeff3 * corr;
	rtd4 += coeff4 * corr;
}

void Particle::position_verlet(double dt, double lx, double ly, Vector G)
{
	rtd2 = (_force + _random_force + _dimple_force) * (1 / _m) + G;
	Vector diff = rtd0 - rtd0_old;
	diff.correct_bc(lx, ly);
	rtd0_new = rtd0 + diff + dt * dt * rtd2;
	rtd1 = (rtd0_new - rtd0_old);
	rtd1.correct_bc(lx, ly);
	rtd1 *= (1 / (2 * dt));

	rtd0_old = rtd0;
	rtd0 = rtd0_new;
}

void Particle::velocity_verlet(double dt, Vector G)
{
	double a1 = dt;
	double a2 = dt * dt / 2;
	
	rtd0 += rtd1 * a1 + rtd2 * a2;
	Vector new_rtd2 = (_force+_random_force) * (1 / _m) + G;
	rtd1 += (rtd2 + new_rtd2) * (dt / 2);
	rtd2 = new_rtd2;
}


void Particle::update_collisions()
{
	if (_contact == 0 && _previous_contact != _contact) {
		_collisions += 1;
	}
	_previous_contact = _contact;
}
