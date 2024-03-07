////////////////////////////////////////////////////////////////////////////
//
// \brief Definition of PFParticle object for LArSoft
//
// \author usher@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/PFParticle.h"

#include <iomanip>
#include <ostream>
#include <utility>

namespace recob {

  PFParticle::PFParticle() : fPdgCode(0), fSelf(0), fParent(PFParticle::kPFParticlePrimary) {}

  PFParticle::PFParticle(int pdgCode,
                         size_t self,
                         size_t parent,
                         const std::vector<size_t>& daughters)
    : fPdgCode(pdgCode), fSelf(self), fParent(parent), fDaughters(daughters)
  {}

  PFParticle::PFParticle(int pdgCode, size_t self, size_t parent, std::vector<size_t>&& daughters)
    : fPdgCode(pdgCode), fSelf(self), fParent(parent), fDaughters(std::move(daughters))
  {}

  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator<<(std::ostream& o, const PFParticle& c)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << "PFParticle hypothesis PDG Code " << std::setw(8) << std::right << c.PdgCode()
      << ", is primary = " << c.IsPrimary() << ", # Daughters " << std::setw(5) << std::right
      << c.NumDaughters() << std::endl;

    return o;
  }

  //----------------------------------------------------------------------
  // less than operator - basically sort in order of index into main collection
  bool operator<(const PFParticle& a, const PFParticle& b)
  {
    return a.Self() < b.Self();
  }

} // namespace
