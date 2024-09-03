#include "lardataobj/MCBase/MCHitCollection.h"

#include <algorithm>

sim::MCHitCollection::MCHitCollection(unsigned int const ch) : fChannel{ch} {}

void sim::MCHitCollection::Reset()
{
  *this = MCHitCollection{};
}

void sim::MCHitCollection::push_back(MCHit const& hit)
{
  bool sort = !empty() && hit < *rbegin();

  std::vector<sim::MCHit>::push_back(hit);
  if (sort) std::sort(begin(), end());
}
