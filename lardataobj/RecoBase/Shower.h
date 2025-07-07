////////////////////////////////////////////////////////////////////////////
//
// \brief Definition of shower object for LArSoft
//
// \author brebel@fnal.gov, modifications by andrzej.szelc@yale.edu and yuntse@slac.stanford.edu
//
////////////////////////////////////////////////////////////////////////////
#ifndef SHOWER_H
#define SHOWER_H

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include <iosfwd>
#include <limits> // std::numeric_limits<>

#include "TVector3.h"

namespace recob {

  class Shower {

  public:
    Shower(); ///Default constructor

  private:
    int fID;
    TVector3 fDCosStart;                   ///< Direction cosines at start of shower
    TVector3 fSigmaDCosStart;              ///< Uncertainty on initial direction cosines
    TVector3 fXYZstart;                    ///< Start point of shower
    TVector3 fSigmaXYZstart;               ///< Uncertainty on start point of shower
    std::vector<double> fTotalEnergy;      ///< Calculated Energy per each plane
    std::vector<double> fSigmaTotalEnergy; ///< Calculated Energy per each plane
    std::vector<double> fdEdx;             ///< Calculated dE/dx per each plane
    std::vector<double> fSigmadEdx;        ///< Calculated dE/dx per each plane

    std::vector<double> fTotalMIPEnergy;      ///< Calculated Energy per each plane
    std::vector<double> fSigmaTotalMIPEnergy; ///< Calculated Energy per each plane
    int fBestPlane;
    /**
     * @brief Shower length [cm].
     *
     * The length of a shower, and should be greater than 0.
     * The characteristic depends on shower reconstruction algorithms.
     * For example, in PCAShowerParticleBuildingAlgorithm (larpandora),
     * it is defined as the three standard deviations of the spacepoint
     * distribution along the principal axis.
     *
     */
    double fLength;

    /**
     * @brief Opening angle [rad].
     *
     * The angle is defined in the @f$ [ 0, \pi/2 ] @f$ range.
     * It is defined as the angle of the shower cone.
     * The characteristic depends on shower reconstruction algorithms.
     * For example, in PCAShowerParticleBuildingAlgorithm (larpandora),
     * it is defined as the ratio of the standard deviation of the
     * spacepoint distribution along the principal axis to that along
     * the secondary axis.
     *
     */
    double fOpenAngle;

    /**
     * @brief The magic constant indicating the invalidity of the
     *        shower length variable
     *
     * This internal constant represents invalid shower length.
     * It is how the function has_length() determines the validity
     * of fLength.
     *
     */
    static constexpr double InvalidLength = std::numeric_limits<double>::lowest();

    /**
     * @brief The magic constant indicating the invalidity of the
     *        opening angle variable
     *
     * This internal constant represents invalid opening angle.
     * It is how the function has_open_angle() determines the validity
     * of fOpenAngle.
     *
     */
    static constexpr double InvalidOpeningAngle = std::numeric_limits<double>::lowest();

  public:
    Shower(TVector3 const& dcosVtx,
           TVector3 const& dcosVtxErr,
           TVector3 const& xyz,
           TVector3 const& xyzErr,
           std::vector<double> TotalEnergy,
           std::vector<double> TotalEnergyErr,
           std::vector<double> dEdx,
           std::vector<double> dEdxErr,
           int bestplane,
           int id,
           double length,
           double openAngle);

    /// Legacy constructor (with no length and no opening angle).
    /// @deprecated Use the complete constructor instead!
    Shower(TVector3 const& dcosVtx,
           TVector3 const& dcosVtxErr,
           TVector3 const& xyz,
           TVector3 const& xyzErr,
           std::vector<double> TotalEnergy,
           std::vector<double> TotalEnergyErr,
           std::vector<double> dEdx,
           std::vector<double> dEdxErr,
           int bestplane,
           int id = util::kBogusI)
      : Shower(dcosVtx,
               dcosVtxErr,
               xyz,
               xyzErr,
               TotalEnergy,
               TotalEnergyErr,
               dEdx,
               dEdxErr,
               bestplane,
               id,
               InvalidLength,
               InvalidOpeningAngle)
    {}

    //set methods
    void set_id(const int id) { fID = id; }
    void set_total_energy(const std::vector<double>& q) { fTotalEnergy = q; }
    void set_total_energy_err(const std::vector<double>& q) { fSigmaTotalEnergy = q; }
    void set_total_MIPenergy(const std::vector<double>& q) { fTotalMIPEnergy = q; }
    void set_total_MIPenergy_err(const std::vector<double>& q) { fSigmaTotalMIPEnergy = q; }
    void set_total_best_plane(const int q) { fBestPlane = q; }

    void set_direction(const TVector3& dir) { fDCosStart = dir; }
    void set_direction_err(const TVector3& dir_e) { fSigmaDCosStart = dir_e; }
    void set_start_point(const TVector3& xyz) { fXYZstart = xyz; }
    void set_start_point_err(const TVector3& xyz_e) { fSigmaXYZstart = xyz_e; }
    void set_dedx(const std::vector<double>& q) { fdEdx = q; }
    void set_dedx_err(const std::vector<double>& q) { fSigmadEdx = q; }
    void set_length(const double& l) { fLength = l; }
    void set_open_angle(const double& a) { fOpenAngle = a; }

    int ID() const;

    const TVector3& Direction() const;
    const TVector3& DirectionErr() const;

    const TVector3& ShowerStart() const;
    const TVector3& ShowerStartErr() const;

    const std::vector<double>& Energy() const;
    const std::vector<double>& EnergyErr() const;

    const std::vector<double>& MIPEnergy() const;
    const std::vector<double>& MIPEnergyErr() const;
    int best_plane() const;
    double Length() const;
    double OpenAngle() const;
    const std::vector<double>& dEdx() const;
    const std::vector<double>& dEdxErr() const;

    //
    // being floating point numbers, equality is a risky comparison;
    // we use anything negative to denote that the following items are not valid
    //

    /// Returns whether the shower has a valid opening angle.
    bool has_open_angle() const;

    /// Returns whether the shower has a valid length.
    bool has_length() const;

    friend std::ostream& operator<<(std::ostream& stream, Shower const& a);

    friend bool operator<(const Shower& a, const Shower& b);

    static_assert(InvalidLength < 0.0f, "Invalid length must be negative!");
    static_assert(InvalidOpeningAngle < 0.0f, "Invalid opening angle must be negative!");

  }; // recob::Shower
}

inline int recob::Shower::ID() const
{
  return fID;
}

inline const TVector3& recob::Shower::Direction() const
{
  return fDCosStart;
}
inline const TVector3& recob::Shower::DirectionErr() const
{
  return fSigmaDCosStart;
}

inline const TVector3& recob::Shower::ShowerStart() const
{
  return fXYZstart;
}
inline const TVector3& recob::Shower::ShowerStartErr() const
{
  return fSigmaXYZstart;
}

inline const std::vector<double>& recob::Shower::Energy() const
{
  return fTotalEnergy;
}
inline const std::vector<double>& recob::Shower::EnergyErr() const
{
  return fSigmaTotalEnergy;
}

inline const std::vector<double>& recob::Shower::MIPEnergy() const
{
  return fTotalMIPEnergy;
}
inline const std::vector<double>& recob::Shower::MIPEnergyErr() const
{
  return fSigmaTotalMIPEnergy;
}
inline int recob::Shower::best_plane() const
{
  return fBestPlane;
}
inline double recob::Shower::Length() const
{
  return fLength;
}
inline double recob::Shower::OpenAngle() const
{
  return fOpenAngle;
}
inline const std::vector<double>& recob::Shower::dEdx() const
{
  return fdEdx;
}
inline const std::vector<double>& recob::Shower::dEdxErr() const
{
  return fSigmadEdx;
}

//
// being floating point numbers, equality is a risky comparison;
// we use anything negative to denote that the following items are not valid
//
inline bool recob::Shower::has_open_angle() const
{
  return fOpenAngle >= 0.0;
}
inline bool recob::Shower::has_length() const
{
  return fLength >= 0.0;
}

#endif // SHOWER_H
