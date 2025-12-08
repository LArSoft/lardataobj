#ifndef TRACKINGTYPE_H
#define TRACKINGTYPE_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "Math/GenVector/AxisAngle.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/SMatrix.h"
#include "TVector3.h"

namespace recob {

  namespace tracking {

    /**
     * Type used for coordinates and values in general.
     * Double32_t is type that matches a 64-bits double when in memory, but is converted to a 32-bit float when written to disk.
     * Given the size and resolution of LArTPC detectors, single-precision floats are sufficient to store the results of tracking
     * algorithms, but it's safer to perform calculations in double precision.
     */
    using Coord_t = Double32_t;

    /// Type for representation of position in physical 3D space. See recob::tracking::Coord_t for more details on the actual type used.
    using Point_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<Coord_t>,
                                                 ROOT::Math::GlobalCoordinateSystemTag>;

    /// Type for representation of momenta in 3D space. See recob::tracking::Coord_t for more details on the actual type used.
    using Vector_t = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Coord_t>,
                                                      ROOT::Math::GlobalCoordinateSystemTag>;

    /// Type of trajectory point list.
    using Positions_t = std::vector<Point_t>;

    /// Type of momentum list.
    using Momenta_t = std::vector<Vector_t>;

    /// Type for representation of space rotations.
    using Rotation_t = ROOT::Math::Rotation3D;

    /// A point in the trajectory, with position and momentum.
    struct TrajectoryPoint_t {

      Point_t position;  ///< position in the trajectory [cm].
      Vector_t momentum; ///< momentum at the trajectory [GeV/c].

      /// Default constructor: sets at origin with no momentum.
      TrajectoryPoint_t() = default;

      /// Constructor: assigns position and momentum.
      TrajectoryPoint_t(Point_t pos, Vector_t mom) : position(pos), momentum(mom) {}

      /// Returns the direction of the trajectory (unit vector of the momentum).
      Vector_t direction() const { return momentum.Unit(); }

    }; // TrajectoryPoint_t

    /// @{
    /// SMatrix and SVector
    using SMatrixSym22 =
      ROOT::Math::SMatrix<Double32_t, 2, 2, ROOT::Math::MatRepSym<Double32_t, 2>>;
    using SMatrixSym33 =
      ROOT::Math::SMatrix<Double32_t, 3, 3, ROOT::Math::MatRepSym<Double32_t, 3>>;
    using SMatrixSym55 =
      ROOT::Math::SMatrix<Double32_t, 5, 5, ROOT::Math::MatRepSym<Double32_t, 5>>;
    using SMatrixSym66 =
      ROOT::Math::SMatrix<Double32_t, 6, 6, ROOT::Math::MatRepSym<Double32_t, 6>>;
    using SMatrix65 = ROOT::Math::SMatrix<Double32_t, 6, 5>;
    using SMatrix56 = ROOT::Math::SMatrix<Double32_t, 5, 6>;
    using SMatrix55 = ROOT::Math::SMatrix<Double32_t, 5, 5>;
    using SMatrix66 = ROOT::Math::SMatrix<Double32_t, 6, 6>;
    using SVector6 = ROOT::Math::SVector<Double32_t, 6>;
    using SVector5 = ROOT::Math::SVector<Double32_t, 5>;
    using SVector3 = ROOT::Math::SVector<Double32_t, 3>;
    using SVector2 = ROOT::Math::SVector<Double32_t, 2>;
    /// @}

  }

}

#endif
