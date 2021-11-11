////////////////////////////////////////////////////////////////////////
// \version $Id:
//
// \brief Definition of ParticleID analysis object
//
// \author brebel@fnal.gov, tjyang@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "lardataobj/AnalysisBase/ParticleID.h"

#include <ostream>

namespace anab{

  //----------------------------------------------------------------------
  ParticleID::ParticleID()
  {
  }

  //----------------------------------------------------------------------
  ParticleID::ParticleID(const std::vector<anab::sParticleIDAlgScores> &ParticleIDAlgScores,
                         const geo::PlaneID& planeID)
  {
    fParticleIDAlgScores = ParticleIDAlgScores;
    fPlaneID = planeID;
  }

  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator<< (std::ostream & o, ParticleID const& a)
  {
    for (size_t i=0; i < a.fParticleIDAlgScores.size(); i++){
      o << "\n ParticleIDAlg " << a.fParticleIDAlgScores.at(i).fAlgName
	<< "\n -- Variable type: " << a.fParticleIDAlgScores.at(i).fVariableType
  << "\n -- Track direction: " << a.fParticleIDAlgScores.at(i).fTrackDir
	<< "\n -- Assuming PDG: " << a.fParticleIDAlgScores.at(i).fAssumedPdg
  << "\n -- Number of degrees of freedom: " << a.fParticleIDAlgScores.at(i).fNdf
	<< "\n -- Value: " << a.fParticleIDAlgScores.at(i).fValue
	<< "\n -- Using planeMask: " << a.fParticleIDAlgScores.at(i).fPlaneMask;
    }
    o << "\n -- Plane id: " << a.fPlaneID;
    o <<std::endl;

    return o;
  }

}
