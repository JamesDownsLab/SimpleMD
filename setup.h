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
	double D; // Noise strength N^2s
};

SystemProperties sysProps{
	1e-6, // timestep
	8.2e-3, // lx 
	4e-3, // ly
	0.0, // x_0
	0.0, // y_0
	1e-8 // D
};

struct MaterialProperties {
	double mass;
	double radius;
	double k; // spring constant
	double e; // coefficient of restitution
};

MaterialProperties matProps{
	0.0387e-3, // mass
	2e-3, // radius
	500e3, // k
	0.2 // e
};

void dump(std::ofstream& os,
	double x, double y, double vx, double vy, int type, MaterialProperties& props) {
	os.precision(8);
	os << x << "\t" << y << "\t" << vx << "\t" << vy << "\t"
		<< props.radius << "\t" << props.mass << "\t" << type << "\t"
		<< props.k << "\t" << props.e << "\t"
		<< "\t0\t0"
		<< "\t0\t0" << "\t0\t0" << "\t0\t0\n";
}

void dump_preamble(std::ofstream& os, SystemProperties& p) {
	os.precision(10);
	os << "#timestep: " << p.timestep << "\n"
		<< "#lx: " << p.lx << "\n"
		<< "#ly: " << p.ly << "\n"
		<< "#x_0: " << p.x_0 << "\n"
		<< "#y_0: " << p.y_0 << "\n"
		<< "#D: " << p.D << "\n";
}

void TwoBalls() {
	std::ofstream fout("initial.random");
	dump_preamble(fout, sysProps);
	dump(fout, 0.00208, 0.002, 2.0, 0.0, 0, matProps);
	dump(fout, 0.00612, 0.002, -2.0, 0.0, 0, matProps);
}

void EnergyInputInvarience(double dt, double D) {
	std::ofstream fout("initial.random");
	SystemProperties props{
		dt,
		49e-3,
		49e-3,
		0.0,
		0.0,
		D
	};
	MaterialProperties mprops{
		0.0387e-3,
		2e-3,
		500e3,
		0.2
	};
	const double PI = 3.1415926;
	//std::random_device rd;
	//std::mt19937_64 gen(rd());
	std::mt19937_64 gen(1);
	std::uniform_real_distribution<double> dis(0.0, 2 * PI);
	dump_preamble(fout, props);
	double v0{ 0.1 };
	for (int i{ 0 }; i < 10; i++) {
		for (int j{ 0 }; j < 10; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, mprops);
		}
	}
}

void TimeStepInvarience(double dt) {
	std::ofstream fout("initial.random");
	SystemProperties props{
		dt,
		49e-3,
		49e-3,
		0.0,
		0.0,
		1e-8
	};
	MaterialProperties mprops{
		0.0387e-3,
		2e-3,
		500e3,
		0.2
	};
	const double PI = 3.1415926;
	//std::random_device rd;
	//std::mt19937_64 gen(rd());
	std::mt19937_64 gen(1);
	std::uniform_real_distribution<double> dis(0.0, 2 * PI);
	dump_preamble(fout, props);
	double v0{ 0.1 };
	for (int i{ 0 }; i < 10; i++) {
		for (int j{ 0 }; j < 10; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, mprops);
		}
	}
}

void HaffsLawTest() {
	std::ofstream fout("initial.random");
	SystemProperties props{
		1e-7,
		49e-3,
		49e-3,
		0.0,
		0.0,
		1e-8
	};
	MaterialProperties mprops{
		0.0387e-3,
		2e-3,
		500e3,
		0.9
	};
	const double PI = 3.1415926;
	//std::random_device rd;
	//std::mt19937_64 gen(rd());
	std::mt19937_64 gen(1);
	std::uniform_real_distribution<double> dis(0.0, 2 * PI);
	dump_preamble(fout, props);
	double v0{ 1.0 };
	for (int i{ 0 }; i < 10; i++) {
		for (int j{ 0 }; j < 10; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, mprops);
		}
	}
}

void BigHaffsLawTest() {
	std::ofstream fout("initial.random");
	SystemProperties props{
		1e-7,
		101e-3,
		101e-3,
		0.0,
		0.0,
		1e-8
	};
	MaterialProperties mprops{
		0.0387e-3,
		2e-3,
		500e3,
		0.2
	};
	const double PI = 3.1415926;
	//std::random_device rd;
	//std::mt19937_64 gen(rd());
	std::mt19937_64 gen(1);
	std::uniform_real_distribution<double> dis(0.0, 2 * PI);
	dump_preamble(fout, props);
	double v0{ 1.0 };
	for (int i{ 0 }; i < 20; i++) {
		for (int j{ 0 }; j < 20; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, mprops);
		}
	}
}


void Balls100() {
	std::ofstream fout("initial.random");
	SystemProperties props{
		1e-7, // timestep
		49e-3, // lx
		49e-3, // ly
		0.0, // x_0
		0.0, // y_0 
		1e-8 // D
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
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, matProps);
		}
	}
}

void Balls10000() {
	std::ofstream fout("initial.random");
	SystemProperties props{
		1e-6,
		499e-3,
		49e-3,
		0.0,
		0.0
	};
	const double PI = 3.1415926;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 2 * PI);
	dump_preamble(fout, props);
	double v0{ 1.0 };
	for (int i{ 0 }; i < 100; i++) {
		for (int j{ 0 }; j < 10; j++) {
			double vx = v0 * cos(dis(gen));
			double vy = v0 * sin(dis(gen));
			vx = 0;
			vy = 0;
			double x = (5.0 * (double)i + 3.0) * 1e-3;
			double y = (5.0 * (double)j + 3.0) * 1e-3;
			dump(fout, x, y, vx, vy, 0, matProps);
		}
	}

}