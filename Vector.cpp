#include "Vector.h"

void Vector::correct_bc(double lx, double ly)
{
	bool result = false;
	if (_x > lx / 2) {
		result = true;
		_x -= lx;
	}
	else if (_x < -lx / 2)
	{
		result = true;
		_x += lx;
	}
	if (_y > ly / 2) {
		result = true;
		_y -= ly;
	}
	else if (_y < -ly / 2) {
		result = true;
		_y += ly;
	}
}
