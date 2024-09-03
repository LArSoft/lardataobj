#ifndef MCWIRECOLLECTION_H
#define MCWIRECOLLECTION_H

// C++ includes
#include "lardataobj/MCBase/MCWire.h"

namespace sim {

  class MCWireCollection : public std::vector<sim::MCWire> {
  public:
    MCWireCollection(unsigned int const ch = sim::kINVALID_UINT);

    unsigned int Channel() const { return fChannel; }
    bool operator<(MCWireCollection const& rhs) const { return fChannel < rhs.fChannel; }

    void Reset();
    void push_back(MCWire const& wire);

  private:
    unsigned int fChannel;
  };
}

// Define a pointer comparison
namespace std {
  template <>
  class less<sim::MCWireCollection*> {
  public:
    bool operator()(sim::MCWireCollection const* lhs, sim::MCWireCollection const* rhs) const
    {
      return *lhs < *rhs;
    }
  };
}

#endif
