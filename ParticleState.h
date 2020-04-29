#pragma once

#include <iostream>

enum class ParticleState {
	Free,
	Confined
};

std::ostream& operator<<(std::ostream& os, const ParticleState& state);

std::istream& operator>>(std::istream& is, ParticleState& state);