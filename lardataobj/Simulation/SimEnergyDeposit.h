////////////////////////////////////////////////////////////////////////
/// \file  lardataobj/Simulation/SimEnergyDeposit.h
/// \brief contains information for a single step in the detector simulation
///
/// \authors Hans Wenzel and William Seligman
////////////////////////////////////////////////////////////////////////
#ifndef LARDATAOBJ_SIMULATION_SIMENERGYDEPOSIT_H
#define LARDATAOBJ_SIMULATION_SIMENERGYDEPOSIT_H

// LArSoft includes
// Define the LArSoft standard geometry types and methods.
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C++ includes
#include <vector>

namespace sim {
  /**
   * @brief Energy deposition in the active material.
   *
   * The detector simulation (presently LArG4, which invokes Geant4)
   * propagates particles through the detector in intervals of "steps".
   * In Geant4, a step is normally defined by the smallest of the distance
   * from the current position of the particle to the point where it
   * enters a new volume boundary, the particle undergoes some "interesting"
   * physics event, or the range of the particle due to its energy falls
   * below a given limit.
   *
   * In `LArG4`, an additional limit is applied: We force the steps to be
   * small (typically 1/10th the wire spacing in the planes of the TPC)
   * so we can process the energy deposited by each step into
   * electron clusters.
   *
   * The SimEnergyDeposit class defines what Geant4 truth information for
   * each step is passed to the ionization -> `sim::SimChannel` conversion,
   * and for the optical-photon -> `sim::SimPhoton` conversion.
   *
   * _William Seligman, Nevis Labs, 10/12/2017_
   */
  class SimEnergyDeposit {
  public:
    // Define the types for the private members below.
    using Length_t = float;
    using Point_t = geo::Point_t;

    // Since we're using LArSoft geometry types, the typical way to
    // construct a SimEnergyDeposit might be:
    //   sim::SimEnergyDeposit sed(numPhotons,
    //                             numElectrons,
    //                             stepEnergy,
    //                             { startX, startY, startZ },
    //                             { endX,   endY,   endZ   },
    //                             startTime,
    //                             endTime,
    //                             trackID,
    //                             pdgCode);

    SimEnergyDeposit(int np = 0,
                     //             int nfp = 0,
                     //             int nsp = 0,
                     int ne = 0,
                     double sy = 0,
                     double e = 0.,
                     geo::Point_t start = {0., 0., 0.},
                     geo::Point_t end = {0., 0., 0.},
                     double t0 = 0.,
                     double t1 = 0.,
                     int id = 0,
                     int pdg = 0,
                     int origTrackID = 0)
      : numPhotons(np)
      //      , numFPhotons(nfp)
      //      , numSPhotons(nsp)
      , numElectrons(ne)
      , scintYieldRatio(sy)
      , edep(e)
      , startPos(start)
      , endPos(end)
      , startTime(t0)
      , endTime(t1)
      , trackID(id)
      , pdgCode(pdg)
      , origTrackID(origTrackID)
    {}

    // Note that even if we store a value as float, we return
    // it as double so the user doesn't have to think about
    // precision issues.

    int NumPhotons() const { return numPhotons; }
    int NumFPhotons() const { return round(numPhotons * scintYieldRatio); }
    int NumSPhotons() const { return round(numPhotons * (1.0 - scintYieldRatio)); }
    int NumElectrons() const { return numElectrons; }
    double ScintYieldRatio() const { return scintYieldRatio; }
    double Energy() const { return edep; }
    geo::Point_t Start() const { return startPos; }
    geo::Point_t End() const { return endPos; }
    double Time() const { return (startTime + endTime) / 2.; }
    int TrackID() const { return trackID; }
    int OrigTrackID() const { return origTrackID; }
    void setTrackID(int id) { trackID = id; }
    int PdgCode() const { return pdgCode; }

    // While it's clear how a SimEnergyDeposit will be created by its
    // constructor, it's not clear how users will want to access its
    // data. So give them as many different kinds of accessors as I
    // can think of.
    geo::Length_t StartX() const { return startPos.X(); }
    geo::Length_t StartY() const { return startPos.Y(); }
    geo::Length_t StartZ() const { return startPos.Z(); }
    double StartT() const { return startTime; }
    geo::Length_t EndX() const { return endPos.X(); }
    geo::Length_t EndY() const { return endPos.Y(); }
    geo::Length_t EndZ() const { return endPos.Z(); }
    double EndT() const { return endTime; }

    // Step mid-point.
    geo::Point_t MidPoint() const { return startPos + 0.5 * (endPos - startPos); }
    geo::Length_t MidPointX() const { return (startPos.X() + endPos.X()) / 2.; }
    geo::Length_t MidPointY() const { return (startPos.Y() + endPos.Y()) / 2.; }
    geo::Length_t MidPointZ() const { return (startPos.Z() + endPos.Z()) / 2.; }
    geo::Length_t X() const { return MidPointX(); }
    geo::Length_t Y() const { return MidPointY(); }
    geo::Length_t Z() const { return MidPointZ(); }
    double T() const { return (startTime + endTime) / 2.; }
    double T0() const { return startTime; }
    double T1() const { return endTime; }
    double E() const { return edep; }

    // Step length. (Recall that the difference between two
    // geo::Point_t objects is a geo::Vector_t; we get the length from
    // spherical coordinates.
    geo::Length_t StepLength() const { return (endPos - startPos).R(); }

    // Just in case someone wants to store sim::SimEnergyDeposit
    // objects in a sorted container, define a sort function. Note
    // that the ideal sort order is dependent of the analysis you're
    // trying to perform; for example, if you're dealing with cosmic
    // rays coming along the y-axis, sorting first by z may cause some
    // tasks like insertions to take a very long time.

    bool operator<(const SimEnergyDeposit& rhs) const
    {
      return trackID < rhs.trackID && startTime < rhs.startTime &&
             startPos.Z() < rhs.startPos.Z() && startPos.Y() < rhs.startPos.Y() &&
             startPos.X() < rhs.startPos.X() && edep > rhs.edep; // sort by _decreasing_ energy
    }

  private:
    // While the accessors above return all values in double
    // precision, store whatever we can in single precision to save
    // memory and disk space.

    // There are roughly 7 digits of decimal precision in a float.
    // This will suffice for energy. A float (as opposed to a double)
    // can hold a little more than 7 digits of decimal precision. The
    // smallest position resolution in the simulation is about 0.1mm,
    // or 10^-4m. With seven digits of precision, that means a float
    // can be accurate to up to the range of 10^3m. That's why the
    // definition of our local Point_t (see above) is based on float,
    // while geo::Point_t is based on double.

    // If the above reasoning is wrong, just change the definition of
    // Length_t near the top of this file. Of course, also edit these
    // comments if you do, because you're a good and responsible
    // programmer.

    // For time, it's possible for long-lived particles like neutrons
    // to deposit energy after billions of ns. Chances are time cuts
    // will take care of that, but let's make sure that any overlay studies
    // won't suffer due to lack of precision.

    int numPhotons;        ///< of scintillation photons
                           //    int           numFPhotons;     ///< of fast scintillation photons
                           //    int           numSPhotons;     ///< of slow scintillation photons
    int numElectrons;      ///< of ionization electrons
    float scintYieldRatio; ///< scintillation yield of LAr
    float edep;            ///< energy deposition (MeV)
    geo::Point_t startPos; ///< positions in (cm)
    geo::Point_t endPos;
    double startTime; ///< (ns)
    double endTime;   ///< (ns)
    int trackID;      ///< simulation track id
    int pdgCode;      ///< pdg code of particle to avoid lookup by particle type later
    int
      origTrackID; ///< complementary simulation track id, kept true to G4 even for shower secondaries/tertiaries etc.
  };
  /*
  // Class utility functions.

  // The format of the sim::SimEnergyDeposit output. I'm using a
  // template for the ostream type, since LArSoft may have some
  // special classes for its output streams.
  template <typename Stream>
  Stream& operator<<(Stream&& os, const sim::SimEnergyDeposit& sed)
  {
    // Note that the geo::Point_t type (returned by Start() and End())
    // has an ostream operator defined for it.
    os << "trackID " << sed.TrackID()
       << " pdgCode=" << sed.PdgCode()
       << " start=" << sed.Start()
       << " t0=" << sed.T0()
       << " end=" << sed.End()
       << " t1=" << sed.T1() << " [cm,ns]"
       << " E=" << sed.E() << "[GeV]"
       << " #photons=" << sed.NumPhotons();
    return os;
  }

  // It can be more memory efficient to sort objects by pointers;
  // e.g., if you've got an unsorted
  // std::vector<sim::SimEnergyDeposit>, create a
  // std::set<sim::SimEnergyDeposit*,sim::CompareSED> so you're not
  // duplicating the objects in memory. The following definition
  // covers sorting the pointers.
  bool compareSED(const SimEnergyDeposit* const lhs, const SimEnergyDeposit* const rhs)
  {
    return (*lhs) < (*rhs);
  }
  */
  typedef std::vector<SimEnergyDeposit> SimEnergyDepositCollection;
} // namespace sim
#endif // LARDATAOBJ_SIMULATION_SIMENERGYDEPOSIT_H
