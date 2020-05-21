#pragma once

#include <vector>

double SQRT3 = sqrt(3);

struct GridResult
{
	double lx;
	double ly;
	std::vector<std::pair<double, double>> coords;
};

GridResult make_grid(double L, int Nx, int Ny) {
	double lx = SQRT3 * Nx * L;
	double ly = Ny * L;
	std::vector<std::pair<double, double>> coords;
	for (int i{ 0 }; i < Nx; i++) {
		for (int j{ 0 }; j < Ny; j++) {
			double x = i * SQRT3 * L;
			double y = j * L;
			coords.push_back({ x, y });
			double x2 = x + L * 0.5 * SQRT3;
			double y2 = y + L * 0.5;
			coords.push_back({ x2, y2 });
		}
	}
	GridResult result{ lx, ly, coords };
	return result;
}