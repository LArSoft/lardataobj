////////////////////////////////////////////////////////////////////////////
// \version
//
// \brief Definition of data product to hold ParticleID information
//
// \author brebel@fnal.gov, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_PARTICLEID_H
#define ANAB_PARTICLEID_H

#include <vector>
#include <iosfwd>
#include <iostream>
#include <iomanip>
#include <bitset>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/ParticleID_VariableTypeEnums.h"

namespace anab {

struct sParticleIDAlgScores { ///< determined particle ID
  std::string fAlgName; ///< Algorithm name (to be defined by experiment). Set to "AlgNameNotSet" by default.
  kVariableType fVariableType; ///< Variable type enum: defined in ParticleID_VariableTypeEnums.h. Set to kNotSet by default.
  kTrackDir fTrackDir; ///< Track direction enum: defined in ParticleID_VariableTypeEnums.h. Set to kNoDirection by default.
  int fNdf; ///< Number of degrees of freedom used by algorithm, if applicable. Set to -9999 by default.
  int fAssumedPdg; ///< PDG of particle hypothesis assumed by algorithm, if applicable. Set to 0 by default.
  float fValue; ///< Result of Particle ID algorithm/test
  std::bitset<8> fPlaneMask; ///< Bitset for PlaneID used by algorithm, allowing for multiple planes and up to 8 total planes. Set to all 0s by default. Convention for bitset is that fPlaneMask[0] (i.e. bit 0) represents plane 0, bit 1 represents plane 1, and so on (with plane conventions defined by the experiment).

  sParticleIDAlgScores(){
  fAlgName = "AlgNameNotSet";
  fVariableType = kNotSet;
  fTrackDir = kNoDirection;
  fAssumedPdg = 0;
  fNdf = -9999;
  fValue = -9999.;
  // fPlaneMask will use default constructor: sets all values to 0
  }

};

  class ParticleID{
  private:
    
    std::vector<sParticleIDAlgScores> fParticleIDAlgScores; ///< Vector of structs to hold outputs from generic PID algorithms
    geo::PlaneID fPlaneID;   ///< Plane id.

  public:

    ParticleID();
    ParticleID(const std::vector<anab::sParticleIDAlgScores> &ParticleIDAlgScores, 
               const geo::PlaneID& planeID);

    friend std::ostream& operator << (std::ostream &o, ParticleID const& a);

    const std::vector<anab::sParticleIDAlgScores> & ParticleIDAlgScores() const;

    const geo::PlaneID& PlaneID() const;

  };

}

inline const std::vector<anab::sParticleIDAlgScores> & anab::ParticleID::ParticleIDAlgScores() const { return fParticleIDAlgScores; }

inline const geo::PlaneID& anab::ParticleID::PlaneID() const {return fPlaneID;}

#endif //ANAB_PARTICLEID_H
