////////////////////////////////////////////////////////////////////////////
//
// Definition of SpacePoint class for LArSoft
//
// SpacePoints are 3D objects that contain pointers to Hits from multiple
// wireplanes that have been identified as matching.
//
// msoderbe@syr.edu
//
////////////////////////////////////////////////////////////////////////////

#ifndef SPACEPOINT_H
#define SPACEPOINT_H

#include <iosfwd>

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include <RtypesCore.h>

namespace recob {

  class SpacePoint {

  public:
    using ID_t = int; ///< type of spacepoint ID

    /// Special value for an invalid ID.
    static constexpr ID_t InvalidID = util::kBogusI;

    SpacePoint(); ///Default constructor

  private:
    ID_t fID;              ///< SpacePoint ID
    Double32_t fXYZ[3];    ///< position of SpacePoint in xyz
    Double32_t fErrXYZ[6]; ///< Error matrix (lower triangular).
    Double32_t fChisq;     ///< Chisquare.

  public:
    SpacePoint(Double32_t const* xyz, Double32_t const* err, Double32_t chisq, int id = InvalidID);

    ID_t ID() const;
    const Double32_t* XYZ() const;
    const Double32_t* ErrXYZ() const;
    Double32_t Chisq() const;

    /// Returns the position of the point in world coordinates [cm]
    geo::Point_t position() const;

    /// Returns the error matrix element for two position coordinates
    /// (`0` for _x_, `1` for _y_ and `2` for _z_)
    double covariance(unsigned int i, unsigned int j) const;

    /// Returns the internal index of correlation structure for coordinates `i` and `j`.
    static constexpr std::size_t covIndex(unsigned int i, unsigned int j);

    friend std::ostream& operator<<(std::ostream& o, const SpacePoint& a);
    friend bool operator<(const SpacePoint& a, const SpacePoint& b);
  };

  /// Comparison of a space point with an ID, for sorting and lookup.
  inline bool operator<(SpacePoint const& s, SpacePoint::ID_t id)
  {
    return s.ID() < id;
  }
  /// Comparison of a space point with an ID, for sorting and lookup.
  inline bool operator<(SpacePoint::ID_t id, SpacePoint const& s)
  {
    return id < s.ID();
  }

}

inline recob::SpacePoint::ID_t recob::SpacePoint::ID() const
{
  return fID;
}
inline const Double32_t* recob::SpacePoint::XYZ() const
{
  return fXYZ;
}
inline const Double32_t* recob::SpacePoint::ErrXYZ() const
{
  return fErrXYZ;
}
inline Double32_t recob::SpacePoint::Chisq() const
{
  return fChisq;
}

inline geo::Point_t recob::SpacePoint::position() const
{
  return geo::Point_t{fXYZ[0], fXYZ[1], fXYZ[2]};
}

#endif //SPACEPOINT_H
