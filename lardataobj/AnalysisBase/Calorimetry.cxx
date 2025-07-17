////////////////////////////////////////////////////////////////////////
//
// \brief Definition of Calorimetry analysis object
//
// \author brebel@fnal.gov, tjyang@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "cetlib_except/exception.h"

#include <iomanip>
#include <iostream>

namespace anab {

  //----------------------------------------------------------------------
  Calorimetry::Calorimetry() : fKineticEnergy(0.), fRange(0.)
  {
    fdEdx.clear();
    fdQdx.clear();
    fResidualRange.clear();
    fDeadWireResR.clear();
    fTrkPitch.clear();
    fXYZ.clear();
    fTpIndices.clear();
  }

  //----------------------------------------------------------------------
  Calorimetry::Calorimetry(float KineticEnergy,
                           std::vector<float> const& dEdx,
                           std::vector<float> const& dQdx,
                           std::vector<float> const& resRange,
                           std::vector<float> const& deadwire,
                           float Range,
                           float TrkPitch,
                           geo::PlaneID planeID)
  {

    fKineticEnergy = KineticEnergy;
    fRange = Range;
    for (size_t i = 0; i != dQdx.size(); ++i) {
      fTrkPitch.push_back(TrkPitch);
      fXYZ.push_back({-999., -999., -999.});
    }
    if (dEdx.size() != resRange.size())
      throw cet::exception("anab::Calorimetry") << "dE/dx and residual range vectors "
                                                << "have different sizes, this is a problem.\n";
    fdEdx.resize(dEdx.size());
    fdQdx.resize(dQdx.size());
    fResidualRange.resize(resRange.size());
    for (size_t i = 0; i < dEdx.size(); ++i) {
      fdEdx[i] = dEdx[i];
      fdQdx[i] = dQdx[i];
      fResidualRange[i] = resRange[i];
    }

    fDeadWireResR.resize(deadwire.size());
    for (size_t i = 0; i < deadwire.size(); ++i) {
      fDeadWireResR[i] = deadwire[i];
    }

    fPlaneID = planeID;
  }

  //----------------------------------------------------------------------
  Calorimetry::Calorimetry(float KineticEnergy,
                           std::vector<float> const& dEdx,
                           std::vector<float> const& dQdx,
                           std::vector<float> const& resRange,
                           std::vector<float> const& deadwire,
                           float Range,
                           std::vector<float> const& TrkPitch,
                           geo::PlaneID planeID)
  {

    fPlaneID = planeID;
    fKineticEnergy = KineticEnergy;
    fRange = Range;
    fTrkPitch = TrkPitch;
    if (dEdx.size() != resRange.size())
      throw cet::exception("anab::Calorimetry") << "dE/dx and residual range vectors "
                                                << "have different sizes, this is a problem.\n";
    for (size_t i = 0; i != dQdx.size(); ++i) {
      fXYZ.push_back({-999., -999., -999.});
    }
    fdEdx.resize(dEdx.size());
    fdQdx.resize(dQdx.size());
    fResidualRange.resize(resRange.size());
    for (size_t i = 0; i < dEdx.size(); ++i) {
      fdEdx[i] = dEdx[i];
      fdQdx[i] = dQdx[i];
      fResidualRange[i] = resRange[i];
    }

    fDeadWireResR.resize(deadwire.size());
    for (size_t i = 0; i < deadwire.size(); ++i) {
      fDeadWireResR[i] = deadwire[i];
    }
  }

  //----------------------------------------------------------------------
  Calorimetry::Calorimetry(float KineticEnergy,
                           std::vector<float> const& dEdx,
                           std::vector<float> const& dQdx,
                           std::vector<float> const& resRange,
                           std::vector<float> const& deadwire,
                           float Range,
                           std::vector<float> const& TrkPitch,
                           std::vector<anab::Point_t> const& XYZ,
                           geo::PlaneID planeID)
    : Calorimetry(KineticEnergy,
                  dEdx,
                  dQdx,
                  resRange,
                  deadwire,
                  Range,
                  TrkPitch,
                  XYZ,
                  std::vector<size_t>(),
                  planeID)
  {}
  //----------------------------------------------------------------------
  Calorimetry::Calorimetry(float KineticEnergy,
                           std::vector<float> const& dEdx,
                           std::vector<float> const& dQdx,
                           std::vector<float> const& resRange,
                           std::vector<float> const& deadwire,
                           float Range,
                           std::vector<float> const& TrkPitch,
                           std::vector<anab::Point_t> const& XYZ,
                           std::vector<size_t> const& TpIndices,
                           geo::PlaneID planeID)
  {

    if (dEdx.size() != resRange.size() || dEdx.size() != dQdx.size() ||
        dEdx.size() != TrkPitch.size() || dEdx.size() != XYZ.size() ||
        (TpIndices.size() > 0 && dEdx.size() != TpIndices.size()))
      throw cet::exception("anab::Calorimetry") << "Input vectors "
                                                << "have different sizes, this is a problem.\n";
    fPlaneID = planeID;
    fKineticEnergy = KineticEnergy;
    fRange = Range;
    fTrkPitch = TrkPitch;
    fdEdx = dEdx;
    fdQdx = dQdx;
    fResidualRange = resRange;
    fXYZ = XYZ;
    fTpIndices = TpIndices;
    fDeadWireResR = deadwire;
  }

  //----------------------------------------------------------------------
  Calorimetry::Calorimetry(float KineticEnergy,
                           std::vector<float> const& dEdx,
                           std::vector<float> const& dQdx,
                           std::vector<float> const& resRange,
                           std::vector<float> const& deadwire,
                           float Range,
                           std::vector<float> const& TrkPitch,
                           std::vector<anab::Point_t> const& XYZ,
                           std::vector<size_t> const& TpIndices,
                           geo::PlaneID planeID,
			   std::vector<float> const& Efield,
			   std::vector<float> const& Phi)
  {

    if (dEdx.size() != resRange.size() || dEdx.size() != dQdx.size() ||
        dEdx.size() != TrkPitch.size() || dEdx.size() != XYZ.size() ||
	dEdx.size() != Efield.size() || dEdx.size() != Phi.size() ||
        (TpIndices.size() > 0 && dEdx.size() != TpIndices.size()))
      throw cet::exception("anab::Calorimetry") << "Input vectors "
                                                << "have different sizes, this is a problem.\n";
    fPlaneID = planeID;
    fKineticEnergy = KineticEnergy;
    fRange = Range;
    fTrkPitch = TrkPitch;
    fdEdx = dEdx;
    fdQdx = dQdx;
    fResidualRange = resRange;
    fXYZ = XYZ;
    fTpIndices = TpIndices;
    fDeadWireResR = deadwire;
    fEfield = Efield;
    fPhi = Phi;
  }

  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator<<(std::ostream& o, Calorimetry const& a)
  {
    o << "Kinetic Energy: " << a.fKineticEnergy << "\n Range: " << a.fRange << std::endl;

    for (size_t n = 0; n < a.fdEdx.size(); ++n)
      o << "dE/dx=" << a.fdEdx[n] << " Residual range=" << a.fResidualRange[n]
        << " dQ/dx=" << a.fdQdx[n] << " (x,y,z)=(" << a.fXYZ[n].X() << "," << a.fXYZ[n].Y() << ","
        << a.fXYZ[n].Z() << ")"
        << " pitch=" << a.fTrkPitch[n] << " planeID=(" << a.fPlaneID.Cryostat << ","
        << a.fPlaneID.TPC << "," << a.fPlaneID.Plane << ")" << std::endl;

    return o;
  }

}
