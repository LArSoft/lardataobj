/**
 * @file    lardataobj/RecoBase/Track.h
 * @brief   Provides `recob::Track` data product.
 * @author  Brian Rebel (brebel@fnal.gov)
 * @ingroup DataProductRecoBase
 * @see     lardataobj/RecoBase/Track.cxx
 */

#ifndef LARDATAOBJ_RECOBASE_TRACK_H
#define LARDATAOBJ_RECOBASE_TRACK_H

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"

#include <iosfwd>
#include <stddef.h>

namespace recob {

  /**
   * @brief Track from a non-cascading particle.
   * @ingroup DataProductRecoBase
   *
   * A `recob::Track` consists of a `recob::TrackTrajectory`, plus additional members relevant for a "fitted" track:
   *
   * * fit &chi;&sup2;
   * * number of degrees of freedom
   * * particle ID hypothesis used in the fit (if any)
   * * covariance matrices at start (vertex) and end positions.
   *
   * Please refer to the `recob::TrackTrajectory` documentation for more information about it;
   * for a discussion on the object type for coordinates see recob::tracking::Coord_t.
   *
   * In terms of interface, `recob::Track` extends `recob::TrackTrajectory`, so that methods of the stored `recob::TrackTrajectory` can be called directly from the `recob::Track interface`,
   * e.g.:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * size_t       NumberTrajectoryPoints()         const { return fTraj.NumberTrajectoryPoints(); }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * Two different parameter conventions are used in a `recob::Track`, and functions to convert from one to the other are provided:
   *
   * 1. Trajectory points and momenta (or directions) are in form of 3-vectors, corresponding to a global Cartesian 6D representation
   * 2. Covariance matrices are stored in a Local 5D representation (so that the covariance matrix is invertible),
   *    where the parameters are defined on the plane orthogonal to the track direction at a given track point.
   *    By construction the local parameters of the track itself are (0,0,0,0,1/p).
   * See `lardataobj/RecoBase/TrackingPlane.h` for more information.
   */
  class Track {

  public:
    using Point_t = tracking::Point_t;
    using Vector_t = tracking::Vector_t;
    using Positions_t = tracking::Positions_t;
    using Momenta_t = tracking::Momenta_t;
    using Rotation_t = tracking::Rotation_t;
    using TrajectoryPoint_t = tracking::TrajectoryPoint_t;

    using SMatrixSym55 = tracking::SMatrixSym55;
    using SMatrixSym66 = tracking::SMatrixSym66;
    using SMatrix65 = tracking::SMatrix65;
    using SMatrix56 = tracking::SMatrix56;
    using SVector6 = tracking::SVector6;
    using SVector5 = tracking::SVector5;

    using PointFlags_t = TrackTrajectory::PointFlags_t;
    using Flags_t = TrackTrajectory::Flags_t;

  protected:
    TrackTrajectory fTraj; ///< Stored trajectory data member
    int fPId = 0;          ///< Particle ID hypothesis used in the fit (if any)
    float fChi2 = -1.;     ///< Fit chi2
    int fNdof = 0.;        ///< Number of degrees of freedom of the fit
    SMatrixSym55
      fCovVertex;         ///< Covariance matrix (local 5D representation) at start point (vertex)
    SMatrixSym55 fCovEnd; ///< Covariance matrix (local 5D representation) at end point
    int fID = -1;         ///< track's ID

  public:
    //Default constructor
    Track() = default;

    Track(TrackTrajectory const& Traj,
          int PId,
          float Chi2,
          int Ndof,
          SMatrixSym55 const& CovVertex,
          SMatrixSym55 const& CovEnd,
          int tkID)
      : fTraj(Traj)
      , fPId(PId)
      , fChi2(Chi2)
      , fNdof(Ndof)
      , fCovVertex(CovVertex)
      , fCovEnd(CovEnd)
      , fID(tkID){};

    Track(TrackTrajectory&& Traj,
          int PId,
          float Chi2,
          int Ndof,
          SMatrixSym55&& CovVertex,
          SMatrixSym55&& CovEnd,
          int tkID)
      : fTraj(std::move(Traj))
      , fPId(PId)
      , fChi2(Chi2)
      , fNdof(Ndof)
      , fCovVertex(std::move(CovVertex))
      , fCovEnd(std::move(CovEnd))
      , fID(tkID){};

    Track(Positions_t&& positions,
          Momenta_t&& momenta,
          Flags_t&& flags,
          bool hasMomenta,
          int PId,
          float Chi2,
          int Ndof,
          SMatrixSym55&& CovVertex,
          SMatrixSym55&& CovEnd,
          int tkID)
      : fTraj(std::move(positions), std::move(momenta), std::move(flags), hasMomenta)
      , fPId(PId)
      , fChi2(Chi2)
      , fNdof(Ndof)
      , fCovVertex(std::move(CovVertex))
      , fCovEnd(std::move(CovEnd))
      , fID(tkID){};

    /// Access to the stored recob::TrackTrajectory
    inline const recob::TrackTrajectory& Trajectory() const { return fTraj; }

    //@{
    /// Various functions related to the presence and the number of (valid) points
    inline size_t NumberTrajectoryPoints() const { return fTraj.NumberTrajectoryPoints(); }
    inline size_t NPoints() const { return fTraj.NPoints(); }
    inline size_t FirstPoint() const { return fTraj.FirstPoint(); }
    inline size_t LastPoint() const { return fTraj.LastPoint(); }
    inline size_t FirstValidPoint() const { return fTraj.FirstValidPoint(); }
    inline size_t NextValidPoint(size_t index) const { return fTraj.NextValidPoint(index); }
    inline size_t PreviousValidPoint(size_t index) const { return fTraj.PreviousValidPoint(index); }
    inline size_t LastValidPoint() const { return fTraj.LastValidPoint(); }
    inline bool HasPoint(size_t i) const { return fTraj.HasPoint(i); }
    inline bool HasValidPoint(size_t i) const { return fTraj.HasValidPoint(i); }
    inline unsigned int CountValidPoints() const { return fTraj.CountValidPoints(); }
    //@}

    //@{
    /// Access to i-th TrajectoryPoint or its Flags
    inline TrajectoryPoint_t TrajectoryPoint(size_t i) const { return fTraj.TrajectoryPoint(i); }
    inline PointFlags_t const& FlagsAtPoint(size_t i) const { return fTraj.FlagsAtPoint(i); }
    //@}

    //@{
    /// Access to track position at different points
    inline Point_t const& Start() const { return fTraj.Start(); }
    inline Point_t const& Vertex() const { return fTraj.Vertex(); }
    inline Point_t const& End() const { return fTraj.End(); }
    inline Point_t const& LocationAtPoint(size_t i) const { return fTraj.LocationAtPoint(i); }
    //@}

    //@{
    /// Access to track direction at different points
    inline Vector_t StartDirection() const { return fTraj.StartDirection(); }
    inline Vector_t VertexDirection() const { return fTraj.VertexDirection(); }
    inline Vector_t EndDirection() const { return fTraj.EndDirection(); }
    inline Vector_t DirectionAtPoint(size_t i) const { return fTraj.DirectionAtPoint(i); }
    //@}

    //@{
    /// Access to track momentum at different points.
    /// The user must check that HasMomentum() returns true to ensure the validity of the result of these functions.
    inline bool HasMomentum() const { return fTraj.HasMomentum(); }
    inline double MomentumAtPoint(unsigned int p) const { return fTraj.MomentumAtPoint(p); }
    inline double VertexMomentum() const { return fTraj.VertexMomentum(); }
    inline double StartMomentum() const { return fTraj.StartMomentum(); }
    inline double EndMomentum() const { return fTraj.EndMomentum(); }
    inline Vector_t const& VertexMomentumVector() const { return fTraj.VertexMomentumVector(); }
    inline Vector_t const& StartMomentumVector() const { return fTraj.StartMomentumVector(); }
    inline Vector_t const& EndMomentumVector() const { return fTraj.EndMomentumVector(); }
    inline Vector_t const& MomentumVectorAtPoint(size_t i) const
    {
      return fTraj.MomentumVectorAtPoint(i);
    }
    //@}

    //@{
    /// Access to covariance matrices
    const SMatrixSym55& StartCovariance() const { return fCovVertex; }
    const SMatrixSym55& VertexCovariance() const { return fCovVertex; }
    const SMatrixSym55& EndCovariance() const { return fCovEnd; }
    //@}

    //@{
    /// Access to position, momentum or covariance at the start and end of the track
    inline std::pair<Point_t, Point_t> Extent() const { return fTraj.Extent(); }
    inline std::pair<Vector_t, Vector_t> Direction() const { return fTraj.Direction(); }
    inline std::pair<SMatrixSym55, SMatrixSym55> Covariances() const
    {
      return std::pair<SMatrixSym55, SMatrixSym55>(fCovVertex, fCovEnd);
    }
    //@}

    //@{
    /// Access to various track properties
    inline double Length(size_t p = 0) const { return fTraj.Length(p); }
    inline float Chi2() const { return fChi2; }
    inline float Chi2PerNdof() const { return fNdof > 0 ? fChi2 / float(fNdof) : util::kBogusF; }
    inline int Ndof() const { return fNdof; }
    inline int ParticleId() const { return fPId; }
    //@}

    //@{
    /// Access to spherical or geographical angles at vertex or at any point
    inline double Theta() const { return fTraj.Theta(); }
    inline double Theta(size_t p) const { return fTraj.Theta(p); }
    inline double Phi() const { return fTraj.Phi(); }
    inline double Phi(size_t p) const { return fTraj.Phi(p); }
    inline double ZenithAngle() const { return fTraj.ZenithAngle(); }
    inline double ZenithAngle(size_t p) const { return fTraj.ZenithAngle(p); }
    inline double AzimuthAngle() const { return fTraj.AzimuthAngle(); }
    inline double AzimuthAngle(size_t p) const { return fTraj.AzimuthAngle(p); }
    //@}

    //@{
    // Calculate rotation matrices between global (x,y,z) and local (u,v,w)
    // coordinate systems based on track direction (fDir).
    // The local w-axis points along the track direction.
    inline Rotation_t GlobalToLocalRotationAtPoint(size_t p) const
    {
      return fTraj.GlobalToLocalRotationAtPoint(p);
    }
    inline Rotation_t LocalToGlobalRotationAtPoint(size_t p) const
    {
      return fTraj.LocalToGlobalRotationAtPoint(p);
    }
    //@}

    //@{
    /// Track ID number, needed to relate a track to its possible track parent (e.g. in case of a refit).
    /// Note that art Assns to the same object are not currently supported.
    /// The < operator is based on the track ID.
    inline int ID() const { return fID; }
    friend bool operator<(const Track& a, const Track& b);
    //@}

    //@{
    /// Accessors to track parameters and covariance matrices in Local5D and Global6D coordinates
    SVector5 VertexParametersLocal5D() const;
    SVector5 EndParametersLocal5D() const;
    const SMatrixSym55& VertexCovarianceLocal5D() const { return fCovVertex; }
    const SMatrixSym55& EndCovarianceLocal5D() const { return fCovEnd; }
    SVector6 VertexParametersGlobal6D() const;
    SVector6 EndParametersGlobal6D() const;
    SMatrixSym66 VertexCovarianceGlobal6D() const;
    SMatrixSym66 EndCovarianceGlobal6D() const;
    //@}

    /// @{
    /// @name Templated version of homonymous functions to access to position, direction, momentum information, and covariances.

    /// Start position. Use e.g. as: @code{.cpp} TVector3 start = track.Start<TVector3>(); @endcode.
    template <typename T>
    inline T Start() const
    {
      return fTraj.Start<T>();
    }

    /// Start position. Use e.g. as: @code{.cpp} TVector3 vertex = track.Vertex<TVector3>(); @endcode.
    template <typename T>
    inline T Vertex() const
    {
      return fTraj.Vertex<T>();
    }

    /// End position. Use e.g. as: @code{.cpp} TVector3 end = track.End<TVector3>(); @endcode.
    template <typename T>
    inline T End() const
    {
      return fTraj.End<T>();
    }

    /// Position at point p. Use e.g. as: @code{.cpp} TVector3 pos = track.LocationAtPoint<TVector3>(p); @endcode.
    template <typename T>
    inline T LocationAtPoint(unsigned int p) const
    {
      return fTraj.LocationAtPoint<T>(p);
    }

    /// Start direction. Use e.g. as: @code{.cpp} TVector3 startdir = track.StartDirection<TVector3>(); @endcode.
    template <typename T>
    inline T StartDirection() const
    {
      return fTraj.StartDirection<T>();
    }

    /// Start direction. Use e.g. as: @code{.cpp} TVector3 vertexdir = track.VertexDirection<TVector3>(); @endcode.
    template <typename T>
    inline T VertexDirection() const
    {
      return fTraj.VertexDirection<T>();
    }

    /// End direction. Use e.g. as: @code{.cpp} TVector3 enddir = track.EndDirection<TVector3>(); @endcode.
    template <typename T>
    inline T EndDirection() const
    {
      return fTraj.EndDirection<T>();
    }

    /// Direction at point p. Use e.g. as: @code{.cpp} TVector3 dir = track.DirectionAtPoint<TVector3>(p); @endcode.
    template <typename T>
    inline T DirectionAtPoint(unsigned int p) const
    {
      return fTraj.DirectionAtPoint<T>(p);
    }

    /// Momentum vector at start point. Use e.g. as: @code{.cpp} TVector3 startmom = track.StartMomentumVector<TVector3>(); @endcode.
    template <typename T>
    inline T StartMomentumVector() const
    {
      return fTraj.StartMomentumVector<T>();
    }

    /// Momentum vector at start point. Use e.g. as: @code{.cpp} TVector3 vertexmom = track.VertexMomentumVector<TVector3>(); @endcode.
    template <typename T>
    inline T VertexMomentumVector() const
    {
      return fTraj.VertexMomentumVector<T>();
    }

    /// Momentum vector at end point. Use e.g. as: @code{.cpp} TVector3 endmom = track.EndMomentumVector<TVector3>(); @endcode.
    template <typename T>
    inline T EndMomentumVector() const
    {
      return fTraj.EndMomentumVector<T>();
    }

    /// Momentum vector at point p. Use e.g. as: @code{.cpp} TVector3 mom = track.MomentumVectorAtPoint<TVector3>(p); @endcode.
    template <typename T>
    inline T MomentumVectorAtPoint(unsigned int p) const
    {
      return fTraj.MomentumVectorAtPoint<T>(p);
    }

    /// Covariance matrix at start point. Use e.g. as: @code{.cpp} TMatrixD startcov = track.StartCovariance<TMatrixD>(); @endcode.
    template <typename T>
    inline T StartCovariance() const;

    /// Covariance matrix at start point. Use e.g. as: @code{.cpp} TMatrixD vertexcov = track.VertexCovariance<TMatrixD>(); @endcode.
    template <typename T>
    inline T VertexCovariance() const
    {
      return StartCovariance<T>();
    }

    /// Covariance matrix at end point. Use e.g. as: @code{.cpp} TMatrixD endcov = track.EndCovariance<TMatrixD>(); @endcode.
    template <typename T>
    inline T EndCovariance() const;

    /// Position at start and end points. Use e.g. as: @code{.cpp} TVector3 start, end; std::tie(start, end) = track.Extent<TVector3>(); @endcode.
    template <typename T>
    inline std::pair<T, T> Extent() const
    {
      return fTraj.Extent<T>();
    }

    /// Direction at start and end points. Use e.g. as: @code{.cpp} TVector3 startdir, enddir; std::tie(startdir, enddir) = track.Direction<TVector3>(); @endcode.
    template <typename T>
    inline std::pair<T, T> Direction() const
    {
      return fTraj.Direction<T>();
    }

    /// Returns a rotation matrix that brings trajectory direction along _z_. Use e.g. as: @code{.cpp} TMatrixD rot = track.GlobalToLocalRotationAtPoint<TMatrixD>(p); @endcode.
    template <typename T>
    inline T GlobalToLocalRotationAtPoint(unsigned int p) const
    {
      return fTraj.GlobalToLocalRotationAtPoint<T>(p);
    }

    /// Returns a rotation matrix bringing relative directions to global. Use e.g. as: @code{.cpp} TMatrixD rot = track.LocalToGlobalRotationAtPoint<TMatrixD>(p); @endcode.
    template <typename T>
    inline T LocalToGlobalRotationAtPoint(unsigned int p) const
    {
      return fTraj.LocalToGlobalRotationAtPoint<T>(p);
    }
    /// @}

  protected:
    friend std::ostream& operator<<(std::ostream& stream, Track const& a);
  };
}

template <typename T>
inline T recob::Track::StartCovariance() const
{
  T result = T(5, 5);
  for (unsigned int i = 0; i < 5; i++) {
    for (unsigned int j = 0; j < 5; j++) {
      result(i, j) = fCovVertex.At(i, j);
    }
  }
  return result;
}

template <typename T>
inline T recob::Track::EndCovariance() const
{
  T result = T(5, 5);
  for (unsigned int i = 0; i < 5; i++) {
    for (unsigned int j = 0; j < 5; j++) {
      result(i, j) = fCovEnd.At(i, j);
    }
  }
  return result;
}

#endif // LARDATAOBJ_RECOBASE_TRACK_H
