#include "ParticleState.h"

std::ostream& operator<<(std::ostream& os, const ParticleState& state) {
  switch (state)
  {
  case ParticleState::Free:
    os << "Free";
    break;

  case ParticleState::Confined:
    os << "Confined";
    break;

  default:
    break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, ParticleState& state) {
  int i;
  is >> i;
  if (i == 0) {
    state = ParticleState::Free;
  }
  else state = ParticleState::Confined;
  return is;
}