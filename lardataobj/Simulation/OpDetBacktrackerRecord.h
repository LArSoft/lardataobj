/*
 * @file   OpDetBacktrackerRecord.h
 * @brief  object containing MC truth information necessary for backtracking Photons.
 * Based on SimChannels by seligman@nevis.columbia.edu
 * @author jason.stock@mines.sdsmt.edu
 * @see    OpDetBacktrackerRecord.cxx
 *
 * This class uses only LArSoft libraries that are header only.
 */

#ifndef LARSIMOBJ_SIMULATION_OPDETBACKTRACKERRECORD_H
#define LARSIMOBJ_SIMULATION_OPDETBACKTRACKERRECORD_H

// C/C++ standard libraries
#include <map>
#include <string>
#include <utility> // std::pair
#include <vector>
namespace sim {

  /// Ionization photons from a Geant4 track
  struct TrackSDP {
    int trackID;      ///< Geant4 supplied trackID
    float energyFrac; ///< fraction of OpHit energy from the particle with this trackID
    float energy;     ///< energy from the particle with this trackID [MeV]

    TrackSDP() {}

    TrackSDP(int id, float phF, float ph) : trackID(id), energyFrac(phF), energy(ph) {}
  };

  /**
   * This class stores information about the Scintillation Light from the simulation
   * of a small step of a track through the TPC active volume.
   *
   * Scintillation information consists of both energy and number of photons.
   * It is of paramount importance to understand what each field stores:
   *
   * - position: where the scintillation occurred (from Geant4 simulation)
   * - track ID: Geant4 track ID of the ionizing particle
   * - energy: amount of energy in G4Step (from Geant4 simulation)
   * - photons: amount of photons reaching the optical detector
   *
   * Note the different definition of the photons respect to the rest:
   * it describes the photons at the OpticlaDetector,
   * while all the other quantities can be related to the moment the scintillation
   * happened.
   *
   * The number of photons is the result of the OpFastScintillation sim, so no real photons are simulated.
   *
   * Each SDP is also typically associated with a time (timePDclock) count, that is
   * the time at which the scintillation photons are generated, as drift time is assumed to be zero for the
   * photons. This time will differ from the time reflected in recob::OpHits as the time stored does not
   * account for any signal processing time. The backtracker will handle that with a FCL parameter
   * called Delay. Each experiment will have to determine the correct Delay to use with their system.
   * For dune10KT the delay is ~300 ns.
   *
   * The SDP (Scintillation Deposited Photons) records are generated in OpFastScintillation.cxx.
   *
   */
  struct SDP {

    /// Type of track ID (the value comes from Geant4)
    typedef int TrackID_t;

    /// Default constructor (sets "bogus" values)
    SDP();

    /// Constructor: copies an SDP, and applies the specified offset to track ID
    SDP(SDP const& sdp, int offset);

    /// Constructor: sets all data members
    SDP(TrackID_t tid, float nPh, float e, float xpos, float ypos, float zpos)
      : trackID(tid), numPhotons(nPh), energy(e), x(xpos), y(ypos), z(zpos)
    {}

    TrackID_t trackID; ///< Geant4 supplied track ID
    float numPhotons;  ///< number of photons at the optical detector for this track ID and time
    float energy;      ///< energy deposited by ionization
    float x;           ///< x position of ionization [cm]
    float y;           ///< y position of ionization [cm]
    float z;           ///< z position of ionization [cm]
  };                   // struct SDP

  /// List of energy deposits at the same time (on this Optical Detector)
  typedef std::pair<int, std::vector<sim::SDP>> timePDclockSDP_t;

  /**
   * @brief Energy deposited on a readout Optical Detector by simulated tracks
   *
   * This class stores a time organized list of scintillation photons detected connected to the G4 tracks
   * they originated from. This class also tracks the energy deposited by those tracks, and what fraction of the
   * energy deposited is realted to the photons detected by this OpDet.
   *
   * The information is organized by time: it is divided by timePDclock ticks (units of ns), and
   * each timePDclock tick where some energy was deposited appears in a separate entry,
   * while the quiet timePDclock ticks are omitted.
   * For each timePDclock, the information is stored as a list of energy deposits;
   * each deposit comes from a single Geant4 track and stores the location where
   * the ionization happened according to the simulation (see `sim::SDP` class).
   *
   * Note that there can be multiple energy deposit records (that is `sim::SDP`)
   * for a single track in a single timePDclock tick.
   */
  class OpDetBacktrackerRecord {
  public:
    /// Type for timePDclock tick used in the internal representation
    typedef timePDclockSDP_t::first_type storedTimePDclock_t;

    /// Type of list of energy deposits for each timePDclock with signal
    // typedef std::vector<timePDclockSDP_t> timePDclockSDPs_t;
    typedef std::map<storedTimePDclock_t, std::vector<sim::SDP>> timePDclockSDPs_t;

  private:
    int iOpDetNum;                     ///< OpticalDetector where the photons were detected
    timePDclockSDPs_t timePDclockSDPs; ///< list of energy deposits for each timePDclock with signal

  public:
    // Default constructor
    OpDetBacktrackerRecord();

    /// Type for iTimePDclock tick used in the interface
    //typedef unsigned short timePDclock_t;
    typedef double timePDclock_t; //This is the G4Time from OpFastScintillation. (ns)

    /// Type of track ID (the value comes from Geant4)
    typedef SDP::TrackID_t TrackID_t;

    /// Constructor: immediately sets the Optical Detector number
    explicit OpDetBacktrackerRecord(int detNum);

    /**
     * @brief Add scintillation photons and energy to this OpticalDetector
     * @param trackID ID of simulated track depositing this energy (from Geant4)
     * @param timePDclock tick when this deposit was collected (ns)
     * @param numberPhotons detected at the OpticalDetector at this time from this track
     * @param xyz coordinates of original location of ionization/scintillation (3D array) [cm]
     * @param energy energy deposited at this point by this track [MeV]
     *
     */
    void AddScintillationPhotons(TrackID_t trackID,
                                 timePDclock_t timePDclock,
                                 double numberPhotons,
                                 double const* xyz,
                                 double energy);

    /// Returns the readout Optical Detector this object describes
    int OpDetNum() const;

    /**
     * @brief Return all the recorded energy deposition within a time interval
     * @param startTimePDclock iTimePDclock tick opening the time window
     * @param endTimePDclock iTimePDclock tick closing the time window (included in the interval)
     * @return a collection of energy deposit information from all tracks
     *
     * This method returns the energy deposited on this Optical Detector by each track
     * ID active in the specified iTimePDclock time interval.
     *
     * Each entry pertains a single track ID. For each entry, all energy
     * deposit information is merged into a single record. It includes:
     * * energy and number of photons, as the integral in the time interval
     * * position, as average weighted by the number of photons
     * * the ID of the track depositing this energy
     *
     * Entries are sorted by track ID number.
     */
    std::vector<sim::SDP> TrackIDsAndEnergies(timePDclock_t startTimePDclock,
                                              timePDclock_t endTimePDclock) const;

    /**
     * @brief Returns all the deposited energy information as stored
     * @return all the deposited energy information as stored in the object
     *
     * The returned list is organized in pairs. Each pair contains all
     * ionization information in a single iTimePDclock tick (collection of `sim::SDP`),
     * and the number of that tick. The information is sorted by increasing timePDclock
     * tick.
     *
     * See the class description for the details of the ionization information
     * content.
     */
    timePDclockSDPs_t const& timePDclockSDPsMap() const;

    /// Returns the total number of scintillation photons on this Optical Detector in the specified timePDclock
    double Photons(timePDclock_t iTimePDclock) const;

    /// Returns the total energy on this Optical Detector in the specified iTimePDclock [MeV]
    double Energy(timePDclock_t iTimePDclock) const;

    // /// Sorts the timeClockSDPs by increasing timePDclock tick
    // void SortTimePDclockSDPs();

    /**
     * @brief Returns energies collected for each track within a time interval
     * @param startTimePDclock iTimePDclock tick opening the time window
     * @param endTimePDclock iTimePDclock tick closing the time window (included in the interval)
     * @return a collection of energy and fraction from each track in interval
     * @see TrackIDsAndEnergies()
     *
     * This method returns the energy deposited on this Optical Detector by each track
     * ID active in the specified iTimePDclock time interval.
     *
     * Each entry pertains a single track ID. For each entry, all energy
     * deposit information is merged into a single record. It includes:
     * * energy of the track, as the integral in the time interval [MeV]
     * * energy fraction respect to the total (see below)
     * * the ID of the track depositing this energy
     *
     * The energy fraction is the energy deposited by the track on this Optical Detector
     * in the specified time interval, divided by the total of the energy
     * deposited by all tracks on this Optical Detector in that same time interval.
     *
     * Entries are sorted by track ID number.
     */
    std::vector<sim::TrackSDP> TrackSDPs(timePDclock_t startTimePDclock,
                                         timePDclock_t endTimePDclock) const;

    /// Comparison: sorts by Optical Detector ID
    bool operator<(const OpDetBacktrackerRecord& other) const;

    /// Comparison: true if OpDetBacktrackerRecords have the same Optical Detector ID
    bool operator==(const OpDetBacktrackerRecord& other) const;

    /**
     * @brief Merges the deposits from another Optical Detector into this one
     * @param opDetNum the sim::OpDetBacktrackerRecord holding information to be merged
     * @param offset track ID offset for the merge
     * @return range of the IDs of the added tracks
     *
     * The information from the specified simulated opDetRecord is added to the
     * current one.
     * This is achieved by appending the energy deposit information (`sim::SDP`)
     * at each iTimePDclock tick from the merged opDetRecord to the list of existing energy
     * deposits for that iTimePDclock tick.
     *
     * In addition, the track IDs of the merged opDetRecord are added an offset,
     * so that they can be distinguished from the existing ones.
     * This is useful when simulating tracks with multiple Geant4 runs. Geant4
     * will reuse track IDs on each run, and using the highest number of track
     * ID from one run as the offset for the next avoids track ID collisions.
     * Note however that this function does not perform any collision check, and
     * it is caller's duty to ensure that the offset is properly large.
     * The return value is a pair including the lowest and the largest track IDs
     * added to this opDetRecord, equivalent to the lowest and the highest track IDs
     * present in the merged opDetRecord, both augmented by the applied offset.
     *
     * The opDetNum number of the merged opDetRecord is ignored.
     */
    std::pair<TrackID_t, TrackID_t> MergeOpDetBacktrackerRecord(
      const OpDetBacktrackerRecord& opDetNum,
      int offset);

    /**
     * @brief Dumps the full content of the OpDetBacktrackerRecord into a stream
     * @tparam Stream an ostream-like stream object
     * @param out the stream to send the information into
     * @param indent indentation of the lines (default: none)
     * @param first_indent indentation for the first line (default: as indent)
     */
    template <typename Stream>
    void Dump(Stream&& out, std::string indent, std::string first_indent) const;

    /// Documentation at `Dump(Stream&&, std::string, std::string) const`.
    template <typename Stream>
    void Dump(Stream&& out, std::string indent = "") const
    {
      Dump(std::forward<Stream>(out), indent, indent);
    }
  };

} // namespace sim

inline bool sim::OpDetBacktrackerRecord::operator<(const sim::OpDetBacktrackerRecord& other) const
{
  return iOpDetNum < other.OpDetNum();
}
inline bool sim::OpDetBacktrackerRecord::operator==(const sim::OpDetBacktrackerRecord& other) const
{
  return iOpDetNum == other.OpDetNum();
}
inline sim::OpDetBacktrackerRecord::timePDclockSDPs_t const&
sim::OpDetBacktrackerRecord::timePDclockSDPsMap() const
{
  return timePDclockSDPs;
}
inline int sim::OpDetBacktrackerRecord::OpDetNum() const
{
  return iOpDetNum;
}

// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class Stream>
void sim::OpDetBacktrackerRecord::Dump(Stream&& out,
                                       std::string indent,
                                       std::string first_indent) const
{
  out << first_indent << "OpDet #" << OpDetNum() << " read " << timePDclockSDPs.size()
      << " timePDclocks:\n";
  double opDet_energy = 0., opDet_photons = 0.;
  for (const auto& timePDclockinfo : timePDclockSDPs) {
    auto const iTimePDclock = timePDclockinfo.first;
    out << indent << "  timePDclock #" << iTimePDclock << " with " << timePDclockinfo.second.size()
        << " SDPs\n";
    double timePDclock_energy = 0., timePDclock_photons = 0.;
    for (const sim::SDP& sdp : timePDclockinfo.second) {
      out << indent << "    (" << sdp.x << ", " << sdp.y << ", " << sdp.z << ") " << sdp.numPhotons
          << " photons, " << sdp.energy << "MeV  (trkID=" << sdp.trackID << ")\n";
      timePDclock_energy += sdp.energy;
      timePDclock_photons += sdp.numPhotons;
    } // for SDPs
    out << indent << "    => timePDclock #" << iTimePDclock << " CH #" << OpDetNum()
        << " collected " << timePDclock_energy << " MeV and " << timePDclock_photons
        << " photons. \n";
    opDet_energy += timePDclock_energy;
    opDet_photons += timePDclock_photons;
  } // for timePDclocks
  out << indent << "  => channel #" << OpDetNum() << " collected " << opDet_photons
      << " photons and " << opDet_energy << " MeV.\n";
} // sim::OpDetBacktrackerRecord::Dump<>()

#endif // LARSIMOBJ_SIMULATION_OPDETBACKTRACKERRECORD_H

////////////////////////////////////////////////////////////////////////
