////////////////////////////////////////////////////////////////////////////
// \version
//
// \brief Definition of data product to hold Calorimetry information
//
// \author brebel@fnal.gov, tjyang@fnal.gov
//
////////////////////////////////////////////////////////////////////////////
#ifndef ANAB_CALORIMETRY_H
#define ANAB_CALORIMETRY_H

#include <iosfwd>
#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

namespace anab {

  using Point_t = recob::tracking::Point_t;

  class Calorimetry {
  public:
    Calorimetry();

    float fKineticEnergy;              ///< determined kinetic energy
    std::vector<float> fdEdx;          ///< dE/dx, should be same size as fResidualRange
    std::vector<float> fdQdx;          ///< dQ/dx
    std::vector<float> fResidualRange; ///< range from end of track
    std::vector<float> fDeadWireResR;  ///< dead wire residual range, collection plane
    float fRange;                      ///< total range of track
    std::vector<float> fTrkPitch;      ///< track pitch on collection plane
    std::vector<Point_t>
      fXYZ; ///< coordinates of space points; for a discussion on the object type for coordinates see recob::tracking::Coord_t.
    std::vector<size_t> fTpIndices; ///< indices of original trajectory points on track
    std::vector<float> fEfield;     ///< E-field strength used for calorimetry
    std::vector<float> fPhi;        ///< track angle at each hit w.r.t. E-field, in degree

  private:
    geo::PlaneID fPlaneID;

  public:
    Calorimetry(float KinematicEnergy,
                std::vector<float> const& dEdx,
                std::vector<float> const& dQdx,
                std::vector<float> const& resRange,
                std::vector<float> const& deadwire,
                float Range,
                float TrkPitch,
                geo::PlaneID planeID);

    Calorimetry(float KineticEnergy,
                std::vector<float> const& dEdx,
                std::vector<float> const& dQdx,
                std::vector<float> const& resRange,
                std::vector<float> const& deadwire,
                float Range,
                std::vector<float> const& TrkPitch,
                geo::PlaneID planeID);

    Calorimetry(float KineticEnergy,
                std::vector<float> const& dEdx,
                std::vector<float> const& dQdx,
                std::vector<float> const& resRange,
                std::vector<float> const& deadwire,
                float Range,
                std::vector<float> const& TrkPitch,
                std::vector<Point_t> const& XYZ,
                geo::PlaneID planeID);

    Calorimetry(float KineticEnergy,
                std::vector<float> const& dEdx,
                std::vector<float> const& dQdx,
                std::vector<float> const& resRange,
                std::vector<float> const& deadwire,
                float Range,
                std::vector<float> const& TrkPitch,
                std::vector<Point_t> const& XYZ,
                std::vector<size_t> const& TpIndices,
                geo::PlaneID planeID);

    Calorimetry(float KineticEnergy,
                std::vector<float> const& dEdx,
                std::vector<float> const& dQdx,
                std::vector<float> const& resRange,
                std::vector<float> const& deadwire,
                float Range,
                std::vector<float> const& TrkPitch,
                std::vector<Point_t> const& XYZ,
                std::vector<size_t> const& TpIndices,
                geo::PlaneID planeID,
                std::vector<float> const& Efield,
                std::vector<float> const& Phi);

    friend std::ostream& operator<<(std::ostream& o, Calorimetry const& a);

    const std::vector<float>& dEdx() const;
    const std::vector<float>& dQdx() const;
    const std::vector<float>& ResidualRange() const;
    const std::vector<float>& DeadWireResRC() const;
    const float& KineticEnergy() const;
    const float& Range() const;
    float TrkPitchC() const;
    const std::vector<float>& TrkPitchVec() const;
    const std::vector<Point_t>& XYZ() const;
    const std::vector<size_t>& TpIndices() const;
    const geo::PlaneID& PlaneID() const;
    const std::vector<float>& Efield() const;
    const std::vector<float>& Phi() const;
  };

}

inline const std::vector<float>& anab::Calorimetry::dEdx() const
{
  return fdEdx;
}
inline const std::vector<float>& anab::Calorimetry::dQdx() const
{
  return fdQdx;
}
inline const std::vector<float>& anab::Calorimetry::ResidualRange() const
{
  return fResidualRange;
}
inline const std::vector<float>& anab::Calorimetry::DeadWireResRC() const
{
  return fDeadWireResR;
}
inline const float& anab::Calorimetry::KineticEnergy() const
{
  return fKineticEnergy;
}
inline const float& anab::Calorimetry::Range() const
{
  return fRange;
}
inline const std::vector<float>& anab::Calorimetry::TrkPitchVec() const
{
  return fTrkPitch;
}
inline float anab::Calorimetry::TrkPitchC() const
{
  if (fTrkPitch.size())
    return fTrkPitch[0];
  else
    return 0;
}
inline const std::vector<anab::Point_t>& anab::Calorimetry::XYZ() const
{
  return fXYZ;
}
inline const std::vector<size_t>& anab::Calorimetry::TpIndices() const
{
  return fTpIndices;
}
inline const geo::PlaneID& anab::Calorimetry::PlaneID() const
{
  return fPlaneID;
}
inline const std::vector<float>& anab::Calorimetry::Efield() const
{
  return fEfield;
}
inline const std::vector<float>& anab::Calorimetry::Phi() const
{
  return fPhi;
}

#endif //ANAB_CALORIMETRY_H
