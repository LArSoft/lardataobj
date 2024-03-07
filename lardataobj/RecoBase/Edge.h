/**
 * @file   Edge.h
 * @brief  An object to define a "edge" which is used to connect
 *         space points in a triangulation algorithm.
 * @author Tracy Usher (usher@slac.stanford.edu)
 * @see    Edge.cxx
 *
 */

#ifndef LARDATAOBJ_RECOBASE_EDGE_H
#define LARDATAOBJ_RECOBASE_EDGE_H

// LArSoft libraries
#include "lardataobj/RecoBase/SpacePoint.h"

// C/C++ libraries
#include <iosfwd> // std::ostream
#include <limits> // std::numeric_limits<>

namespace recob {

  /**
     * @brief Edge is an object containing the results of a Principal Components
     *        Analysis of a group of space points.
     *
     * The edge contains references to an emanating space point
     * (`FirstPointID()`) and to an ending one (`SecondPointID()`).
     * For convenience, it also stores the distance between those points.
     *
     * To look up for a referenced space point, the easiest way is to start
     * from a sorted list of space points (`recob::SpacePoint` sorts by ID
     * value):
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * if (!std::is_sorted(points.begin(), points.end()))
     *   throw std::runtime_error("Space points not sorted!");
     *
     * // find the first space point
     * auto const iFirstPoint = std::lower_bound
     *   (points.begin(), points.end(), edge.FirstPointID());
     *
     * if ((iFirstPoint == points.end()) || (iFirstPoint->ID() != edge.FirstPointID())) {
     *   throw std::runtime_error
     *     ("First point not found: ID=" + std::to_string(edge.FirstPointID()));
     * }
     * recob::SpacePoint const& firstPoint = *iFirstPoint;
     *
     * // find the second space point
     * auto const iSecondPoint = std::lower_bound
     *   (points.begin(), points.end(), edge.SecondPointID());
     *
     * if ((iSecondPoint == points.end()) || (iSecondPoint->ID() != edge.SecondPointID())) {
     *   throw std::runtime_error
     *     ("Second point not found: ID=" + std::to_string(edge.SecondPointID()));
     * }
     * recob::SpacePoint const& secondPoint = *iSecondPoint;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     */
  class Edge {
  public:
    /// Type to represent `recob::Edge` IDs.
    using ID_t = unsigned int;

    /// Type to represent `recob::SpacePoint` IDs.
    using SpacePointID_t = recob::SpacePoint::ID_t;

    /// Special value for an invalid edge ID.
    static constexpr ID_t InvalidID = std::numeric_limits<ID_t>::max();

    /// Default constructor (all invalid IDs).
    Edge() = default;

    /**
         * @brief Constructor: assigns all values.
         * @param length the length of the edge [cm]
         * @param firstPointID ID of the emanating space point
         * @param secondPointID ID of the ending space point
         * @param id _(default: `InvalidID`) ID of this edge
         *
         */
    Edge(const double length,
         SpacePointID_t firstPointID,
         SpacePointID_t secondPointID,
         ID_t id = InvalidID);

    /**
         * @brief Constructor: uses the specified spacepoints.
         * @param firstPoint the emanating space point
         * @param secondPoint the ending space point
         * @param id _(default: `InvalidID`) ID of this edge
         *
         */
    Edge(SpacePoint const& firstPoint, SpacePoint const& secondPoint, ID_t id = InvalidID);

    /// @{
    /// @name Access

    /// Returns the length of this edge [cm].
    double Length() const { return fLength; }

    /// Returns the ID of the SpacePoint this edge emanates from.
    SpacePointID_t FirstPointID() const { return fFirstPointID; }

    /// Returns the ID of the SpacePoint this edge ends on.
    SpacePointID_t SecondPointID() const { return fSecondPointID; }

    /// Returns the ID of this edge.
    ID_t ID() const { return fID; }

    /// @}

  private:
    /// Length of this Edge [cm].
    double fLength = 0.0;
    /// ID of the SpacePoint edge emanates from.
    SpacePointID_t fFirstPointID = recob::SpacePoint::InvalidID;
    /// ID of the SpacePoint edge ends on.
    SpacePointID_t fSecondPointID = recob::SpacePoint::InvalidID;
    /// Edge ID.
    ID_t fID = InvalidID;

  }; // class Edge

  /// Streaming operator: prints the edge `a` into the specified stream.
  std::ostream& operator<<(std::ostream& o, const Edge& a);

  /// Comparison operator: strict ordering of edge by ID.
  inline bool operator<(const Edge& a, const Edge& b)
  {
    return a.ID() < b.ID();
  }
  /// Comparison operator: strict ordering of edge by ID.
  inline bool operator<(const Edge& e, Edge::ID_t id)
  {
    return e.ID() < id;
  }
  /// Comparison operator: strict ordering of edge by ID.
  inline bool operator<(Edge::ID_t id, const Edge& e)
  {
    return id < e.ID();
  }

} // namespace recob

#endif // LARDATAOBJ_RECOBASE_EDGE_H
