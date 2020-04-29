#pragma once
#include <fstream>
#include <random>
#include <math.h>

struct SystemProperties {
	double timestep;
	double lx;
	double ly;
	double x_0;
	double y_0;
};

SystemProperties sysProps{
	1e-7,
	8.2e-3,
	4e-3,
	0.0,
	0.0
};

struct MaterialProperties {
	double mass;
	double radius;
	double youngs;
	double mu;
	double A;
	double gamma;
};

MaterialProperties matProps{
	0.0387e-3,
	2e-3,
	7.32e6,
	0.5,
	3e-5,
	2.06
};

void dump(std::ofstream& os,
	double x, double y, double vx, double vy, int type, MaterialProperties& props) {
	os.precision(5);
	os << x << "\t" << y << "\t" << vx << "\t" << vy << "\t"
		<< props.radius << "\t" << props.mass << "\t" << type << "\t"
		<< props.youngs << "\t" << props.A << "\t" << props.mu << "\t" << props.gamma
		<< "\t0\t0"
		<< "\t0\t0" << "\t0\t0" << "\t0\t0\n";
}

void dump_preamble(std::ofstream& os, SystemProperties& p) {
	os << "#timestep: " << p.timestep << "\n"
		<< "#lx: " << p.lx << "\n"
		<< "#ly: " << p.ly << "\n"
		<< "#x_0: " << p.x_0 << "\n"
		<< "#y_0: " << p.y_0 << "\n";
}

void TwoBalls() {
	std::ofstream fout("initial.random");
	dump_preamble(fout, sysProps);
	dump(fout, 0.00208, 0.002, 2.0, 0.0, 0, matProps);
	dump(fout, 0.00612, 0.002, -2.0, 0.0, 0, matProps);
}

void Balls100() {
	std::ofstream fout("initial.random");
	SystemProperties props{
		1e-6,
		49e-3,
		49e-3,
		0.0,
		0.0
	};
	const double PI = 3.1415926;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 2*PI);
	dump_preamble(fout, props);
	double v0{ 1.0 };
	for (int i{ 0 }; i < 10; i++) {
		for (int j{ 0 }; j < 10; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			vx = 0;
			vy = 0;
			dump(fout, (5 * i + 3) * 1e-3, (5 * j + 3) * 1e-3, vx, vy, 0, matProps);
		}
	}

}