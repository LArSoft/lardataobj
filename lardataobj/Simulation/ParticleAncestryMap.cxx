#ifndef PARTICLEANCESTRYMAP_CXX
#define PARTICLEANCESTRYMAP_CXX

#include "lardataobj/Simulation/ParticleAncestryMap.h"

#include <limits>

namespace sim {

  void ParticleAncestryMap::SetMap(const std::map<int, std::set<int>>& map)
  {
    fParticleMap = map;
  }

  std::map<int, std::set<int>> const& ParticleAncestryMap::GetMap() const
  {
    return fParticleMap;
  }

  bool ParticleAncestryMap::HasDroppedDescendants(const int trackid) const
  {
    return fParticleMap.count(trackid) != 0;
  }

  std::set<int> const& ParticleAncestryMap::GetAllDroppedDescendants(const int trackid) const
  {
    return fParticleMap.at(trackid);
  }

  int ParticleAncestryMap::GetAncestor(const int trackid) const
  {
    for (auto const& [ancestor, descendants] : fParticleMap) {
      if (descendants.count(trackid) != 0) return ancestor;
    }

    return -std::numeric_limits<int>::max();
  }

  bool ParticleAncestryMap::Exists(const int trackid) const
  {
    return trackid != -std::numeric_limits<int>::max();
  }
}

#endif
