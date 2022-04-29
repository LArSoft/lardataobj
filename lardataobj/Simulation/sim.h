/// \file lardataobj/Simulation/sim.h
///
/// \brief Tools and modules for checking out the basics of the Monte Carlo
///
/// \author brebel@fnal.gov

#ifndef LARDATAOBJ_SIMULATION_SIM_H
#define LARDATAOBJ_SIMULATION_SIM_H

#include <limits>
#include "TRandom3.h"

///Monte Carlo Simulation
namespace sim{

  unsigned int GetRandomNumberSeed();

  // any track id method returns sim::Particle:NoParticleId, it means the
  // associated particle was too low-energy to be written by the
  // detector Monte Carlo.
  static const int NoParticleId = std::numeric_limits<int>::min();

}

inline unsigned int sim::GetRandomNumberSeed(){

    // the maximum allowed seed for the art::RandomNumberGenerator
    // is 900000000. Use TRandom3 to get the seed value in that range.
    // Instantiating TRandom3 with a 0 means that its seed is set based
    // on the TUUID and should always be random, even for jobs running on the
    // same machine
    TRandom3 rand(0);
    return rand.Integer(900000000);
}


#endif// LARDATAOBJ_SIMULATION_SIM_H
