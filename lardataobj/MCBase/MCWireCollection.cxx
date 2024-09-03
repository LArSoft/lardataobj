#include "lardataobj/MCBase/MCWireCollection.h"

#include <algorithm>

sim::MCWireCollection::MCWireCollection(unsigned int const ch) : fChannel{ch} {}

void sim::MCWireCollection::Reset()
{
  *this = MCWireCollection();
}

void sim::MCWireCollection::push_back(MCWire const& wire)
{
  bool sort = !empty() && wire < *rbegin();

  std::vector<sim::MCWire>::push_back(wire);
  if (sort) std::sort(begin(), end());
}
