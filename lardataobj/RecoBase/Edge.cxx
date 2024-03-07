/**
 * @file   Edge.cxx
 * @brief  Definition of Principal Components Axis object for LArSoft.
 * @author Tracy Usher (usher@slac.stanford.edu)
 * @see    Edge.h
 *
 */

#include "lardataobj/RecoBase/Edge.h"

// C/C++ standard libraries
#include <cmath> // std::sqrt()
#include <iomanip>
#include <ostream>

namespace {

  double sqr(double v)
  {
    return v * v;
  }

  double spacePointDistance(recob::SpacePoint const& a, recob::SpacePoint const& b)
  {
    double const* aXYZ = a.XYZ();
    double const* bXYZ = b.XYZ();
    return std::sqrt(sqr(aXYZ[0] - bXYZ[0]) + sqr(aXYZ[1] - bXYZ[1]) + sqr(aXYZ[2] - bXYZ[2]));
  } // spacePointDistance()

} // local namespace

namespace recob {

  //----------------------------------------------------------------------
  Edge::Edge(const double length,
             SpacePointID_t firstPointID,
             SpacePointID_t secondPointID,
             ID_t id /* = InvalidID */)
    : fLength(length), fFirstPointID(firstPointID), fSecondPointID(secondPointID), fID(id)
  {}

  //----------------------------------------------------------------------
  Edge::Edge(SpacePoint const& firstPoint, SpacePoint const& secondPoint, ID_t id /* = InvalidID */)
    : fLength(::spacePointDistance(firstPoint, secondPoint))
    , fFirstPointID(firstPoint.ID())
    , fSecondPointID(secondPoint.ID())
    , fID(id)
  {}

  //----------------------------------------------------------------------
  // ostream operator.
  //
  std::ostream& operator<<(std::ostream& o, const Edge& a)
  {
    o << std::setiosflags(std::ios::fixed) << std::setprecision(2);
    o << " Edge ID " << a.ID() << " has length " << std::setw(6) << a.Length() << " cm"
      << std::endl;
    o << "   - First Point ID: " << a.FirstPointID() << ", second point ID: " << a.SecondPointID()
      << std::endl;

    return o;
  }

} // End of namespace
