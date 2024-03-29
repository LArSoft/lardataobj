#ifndef TRACKINGPLANE_H
#define TRACKINGPLANE_H

#include "lardataobj/RecoBase/TrackingTypes.h"

namespace recob {
  namespace tracking {

    /// \file  lardataobj/RecoBase/TrackingPlane.h
    /// \class recob::tracking::Plane
    ///
    /// \brief Class defining a plane for tracking.
    ///        It provides various functionalities to convert track parameters and covariance matrices from global to local coordinates.
    ///
    /// \author  G. Cerati (FNAL, MicroBooNE)
    /// \date    2017
    /// \version 1.0
    ///
    /// The plane is constructed from a point in global coordinates, and a direction unit vector orthogonal to the plane.
    /// For instance, they may refer to the position and direction of a trajectory at a given point,
    /// or to the center of a wire and the direction orthogonal to both the wire and drift directions.
    ///
    /// The plane and parameters are defined along the lines of trkf::SurfXYZPlane, but with a different notation for the angles
    /// which are now called alpha and beta to avoid confusion with theta and phi in global spherical coordinates:
    /// (alpha, beta) correspond to (alpha := Rotation angle about y'-axis (projected Lorentz angle) = atan2(nx, hypot(ny, nz)),
    ///                               beta := Rotation angle about x-axis (wire angle) = atan2(-ny, nz))
    /// The 3D local positions are defined in terms of the unit vectors u,v,w where w is along the normal direction to the plane,
    /// u is perpendicular to w and coplanar to the global x axis and with position sign along the positive sign of x,
    /// v forms a right handed orthonormal basis with u and w.
    /// The 5D track parameters are u,v,du/dw,dv/dw,1/p:
    /// u,v are the local positions, du, dv, and dw are the local directions, and 1/p is the inverse of the track momentum.
    /// The global 6D coordinates are formed by the track position and momentum (or direction) at a given point.
    ///
    /// The discussion above refers to a detector with drift direction along the x axis;
    /// this class will have to be extended for a detector with a different drift direction.

    class Plane {

      ///
      /// Struct caching trigonometric function results.
      struct TrigCache {
      public:
        TrigCache(const Vector_t& planeDir)
        {
          const double diryz = std::hypot(planeDir.Y(), planeDir.Z());
          fCosA = diryz;
          fSinA = planeDir.X();
          fCosB = (diryz != 0.0) ? planeDir.Z() / diryz : 1.0;
          fSinB = (diryz != 0.0) ? -planeDir.Y() / diryz : 0.0;
        }
        double fCosA;
        double fSinA;
        double fCosB;
        double fSinB;
      };

    public:
      ///
      /// Constructor from reference position on the plane and direction orthogonal to the plane
      Plane(const Point_t& planePos, const Vector_t& planeDir)
        : fPlanePos(planePos), fPlaneDir(planeDir.Unit()), fTrigCache(planeDir.Unit())
      {}

      ///
      /// Reference position on the plane
      Point_t const& position() const { return fPlanePos; }

      ///
      /// Reference direction orthogonal to the plane
      Vector_t const& direction() const { return fPlaneDir; }

      //@{
      /// Function to convert parameters from global to local coordinates. Local coordinates are on the plane with origin at planePos and orthogonal to planeDir. It is responsibility of the user to make sure the global position lies on the plane
      inline SVector5 Global6DToLocal5DParameters(const SVector6& par6d) const
      {
        return Global6DToLocal5DParameters(par6d, fPlanePos, fPlaneDir, fTrigCache);
      }
      SVector5 Global6DToLocal5DParameters(const SVector6& par6d,
                                           const Point_t& planePos,
                                           const Vector_t& planeDir,
                                           const TrigCache& trigCache) const;
      //@}

      //@{
      /// Function to convert parameters from local to global coordinates. Local coordinates are on the plane with origin at planePos and orthogonal to planeDir. trackAlongPlaneDir is as given by trackDir.Dot(planeDir)>0
      inline SVector6 Local5DToGlobal6DParameters(const SVector5& par5d,
                                                  bool trackAlongPlaneDir = true) const
      {
        return Local5DToGlobal6DParameters(
          par5d, fPlanePos, fPlaneDir, fTrigCache, trackAlongPlaneDir);
      }
      SVector6 Local5DToGlobal6DParameters(const SVector5& par5d,
                                           const Point_t& planePos,
                                           const Vector_t& planeDir,
                                           const TrigCache& trigCache,
                                           bool trackAlongPlaneDir = true) const;
      //@}

      //@{
      /// Compute the jacobian to translate track covariance from local to global coordinates. The track momentum (or direction) is needed to compute the jacobian. Local coordinates are on the plane orthogonal to planeDir (it may be the same direction as the momentum, but the function is generic).
      inline SMatrix65 Local5DToGlobal6DJacobian(bool hasMomentum,
                                                 const Vector_t& trackMomOrDir) const
      {
        return Local5DToGlobal6DJacobian(hasMomentum, trackMomOrDir, fPlaneDir, fTrigCache);
      }
      SMatrix65 Local5DToGlobal6DJacobian(bool hasMomentum,
                                          const Vector_t& trackMomOrDir,
                                          const Vector_t& planeDir,
                                          const TrigCache& trigCache) const;
      //@}

      //@{
      /// Compute the jacobian to translate track covariance from global to local coordinates. The track momentum (or direction) is needed to compute the jacobian. Local coordinates are on the plane orthogonal to planeDir (it may be the same direction as the momentum, but the function is generic). Warning: some information may be lost in degenerate cases, e.g. the unceratinty along z position when converting to a x-y plane (fixed z)
      inline SMatrix56 Global6DToLocal5DJacobian(bool hasMomentum,
                                                 const Vector_t& trackMomOrDir) const
      {
        return Global6DToLocal5DJacobian(hasMomentum, trackMomOrDir, fPlaneDir, fTrigCache);
      }
      SMatrix56 Global6DToLocal5DJacobian(bool hasMomentum,
                                          const Vector_t& trackMomOrDir,
                                          const Vector_t& planeDir,
                                          const TrigCache& trigCache) const;
      //@}

      /// Translate track covariance from local to global coordinates. The track momentum (or direction) is needed to compute the jacobian. Local coordinates are on the plane orthogonal to planeDir (it may be the same direction as the momentum, but the function is generic).
      inline SMatrixSym66 Local5DToGlobal6DCovariance(SMatrixSym55 cov5d,
                                                      bool hasMomentum,
                                                      const Vector_t& trackMomOrDir) const
      {
        return Local5DToGlobal6DCovariance(cov5d, hasMomentum, trackMomOrDir, fPlaneDir);
      }

      /// Translate track covariance from global to local coordinates. The track momentum (or direction) is needed to compute the jacobian. Local coordinates are on the plane orthogonal to planeDir (it may be the same direction as the momentum, but the function is generic). Warning: some information may be lost in degenerate cases, e.g. the unceratinty along z position when converting to a x-y plane (fixed z)
      inline SMatrixSym55 Global6DToLocal5DCovariance(SMatrixSym66 cov6d,
                                                      bool hasMomentum,
                                                      const Vector_t& trackMomOrDir) const
      {
        return Global6DToLocal5DCovariance(cov6d, hasMomentum, trackMomOrDir, fPlaneDir);
      }

      //@{
      /// Calculate rotation matrices from global (x,y,z) to local (u,v,w) coordinates
      inline Rotation_t Global3DToLocal3DRotation() const
      {
        return Global3DToLocal3DRotation(fPlaneDir, fTrigCache);
      }
      Rotation_t Global3DToLocal3DRotation(const Vector_t& planeDir,
                                           const TrigCache& trigCache) const;
      //@}

      //@{
      /// Calculate rotation matrices from local (u,v,w) to global (x,y,z) coordinates
      inline Rotation_t Local3DToGlobal3DRotation() const
      {
        return Local3DToGlobal3DRotation(fPlaneDir, fTrigCache);
      }
      Rotation_t Local3DToGlobal3DRotation(const Vector_t& planeDir,
                                           const TrigCache& trigCache) const;
      //@}

      //@{
      /// Return cached values of trigonometric function for angles defining the plane
      double cosAlpha() const { return fTrigCache.fCosA; }
      double sinAlpha() const { return fTrigCache.fSinA; }
      double cosBeta() const { return fTrigCache.fCosB; }
      double sinBeta() const { return fTrigCache.fSinB; }
      //@}

      /// \copydoc Plane::Global6DToLocal5DParameters
      static SVector5 Global6DToLocal5DParameters(const SVector6& par6d,
                                                  const Point_t& planePos,
                                                  const Vector_t& planeDir)
      {
        Plane p(planePos, planeDir);
        return p.Global6DToLocal5DParameters(par6d);
      }
      /// \copydoc Plane::Local5DToGlobal6DParameters
      static SVector6 Local5DToGlobal6DParameters(const SVector5& par5d,
                                                  const Point_t& planePos,
                                                  const Vector_t& planeDir,
                                                  bool trackAlongPlaneDir = true)
      {
        Plane p(planePos, planeDir);
        return p.Local5DToGlobal6DParameters(par5d, trackAlongPlaneDir);
      }
      /// \copydoc Plane::Local5DToGlobal6DJacobian
      static SMatrix65 Local5DToGlobal6DJacobian(bool hasMomentum,
                                                 const Vector_t& trackMomOrDir,
                                                 const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Local5DToGlobal6DJacobian(hasMomentum, trackMomOrDir);
      }
      /// \copydoc Plane::Local5DToGlobal6DJacobian
      static SMatrix65 Local5DToGlobal6DJacobian(bool hasMomentum,
                                                 const SVector6& par6d,
                                                 const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Local5DToGlobal6DJacobian(hasMomentum, Vector_t(par6d[3], par6d[4], par6d[5]));
      }
      /// \copydoc Plane::Global6DToLocal5DJacobian
      static SMatrix56 Global6DToLocal5DJacobian(bool hasMomentum,
                                                 const Vector_t& trackMomOrDir,
                                                 const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Global6DToLocal5DJacobian(hasMomentum, trackMomOrDir);
      }
      /// \copydoc Plane::Global6DToLocal5DJacobian
      static SMatrix56 Global6DToLocal5DJacobian(bool hasMomentum,
                                                 const SVector6& par6d,
                                                 const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Global6DToLocal5DJacobian(hasMomentum, Vector_t(par6d[3], par6d[4], par6d[5]));
      }
      /// \copydoc Plane::Local5DToGlobal6DCovariance
      static SMatrixSym66 Local5DToGlobal6DCovariance(SMatrixSym55 cov5d,
                                                      bool hasMomentum,
                                                      const Vector_t& trackMomOrDir,
                                                      const Vector_t& planeDir)
      {
        return ROOT::Math::Similarity(
          Local5DToGlobal6DJacobian(hasMomentum, trackMomOrDir, planeDir), cov5d);
      }
      /// \copydoc Plane::Global6DToLocal5DCovariance
      static SMatrixSym55 Global6DToLocal5DCovariance(SMatrixSym66 cov6d,
                                                      bool hasMomentum,
                                                      const Vector_t& trackMomOrDir,
                                                      const Vector_t& planeDir)
      {
        return ROOT::Math::Similarity(
          Global6DToLocal5DJacobian(hasMomentum, trackMomOrDir, planeDir), cov6d);
      }
      /// \copydoc Plane::Global3DToLocal3DRotation
      static Rotation_t Global3DToLocal3DRotation(const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Global3DToLocal3DRotation();
      }
      /// \copydoc Plane::Local3DToGlobal3DRotation
      static Rotation_t Local3DToGlobal3DRotation(const Vector_t& planeDir)
      {
        Plane p(Point_t(), planeDir);
        return p.Local3DToGlobal3DRotation();
      }

    private:
      Point_t fPlanePos;    ///< Position of a point on the plane
      Vector_t fPlaneDir;   ///< Direction vector othogonal to the plane
      TrigCache fTrigCache; ///< Cached trigonometric function values
    };

  }
}

#endif
