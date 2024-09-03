#ifndef MCHITCOLLECTION_H
#define MCHITCOLLECTION_H

// C++ includes
#include "lardataobj/MCBase/MCHit.h"
#include <vector>

namespace sim {

  class MCHitCollection : public std::vector<sim::MCHit> {
  public:
    MCHitCollection(const unsigned int ch = ::sim::kINVALID_UINT);

    void Reset();
    void push_back(MCHit const& hit);

    unsigned int Channel() const { return fChannel; }
    bool operator<(MCHitCollection const& rhs) const { return fChannel < rhs.fChannel; }

  private:
    unsigned int fChannel;
  };
}

// Define a pointer comparison
namespace std {
  template <>
  class less<sim::MCHitCollection*> {
  public:
    bool operator()(sim::MCHitCollection const* lhs, sim::MCHitCollection const* rhs) const
    {
      return *lhs < *rhs;
    }
  };
}

#endif
