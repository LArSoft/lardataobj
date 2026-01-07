/**
 * @file   SimChannel.h
 * @brief  object containing MC truth information necessary for making
 *         RawDigits and doing back tracking
 * @author seligman@nevis.columbia.edu
 * @see    SimChannel.cxx
 *
 * This class uses only LArSoft libraries that are header only.
 */

#ifndef LARDATAOBJ_SIMULATION_SIMCHANNEL_H
#define LARDATAOBJ_SIMULATION_SIMCHANNEL_H

// LArSoftObj libraries
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::kBogusI
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"          // raw::ChannelID_t

// C/C++ standard libraries
#include <string>
#include <utility> // std::pair
#include <vector>

namespace sim {

  /// Ionization energy from a Geant4 track
  struct TrackIDE {
    int trackID;        ///< Geant4 supplied trackID
    float energyFrac;   ///< fraction of hit energy from the particle with this trackID
    float energy;       ///< energy from the particle with this trackID [MeV]
    float numElectrons; ///< number of electrons from the particle detected on the wires
    int
      origTrackID; ///< Geant4 supplied trackID, including no modification for shower secondaries/tertiaries

    TrackIDE() {}

    TrackIDE(int id, float ef, float e, float ne, int gid = util::kBogusI)
      : trackID(id), energyFrac(ef), energy(e), numElectrons(ne), origTrackID(gid)
    {}
  };

  /**
   * @brief Ionization at a point of the TPC sensitive volume
   *
   * This class stores information about the ionization from the simulation
   * of a small step of a track through the TPC active volume.
   *
   * Ionization information consists of both energy and number of electrons.
   * It is of paramount importance to understand what each field stores:
   *
   * - position: where the ionization occurred (from Geant4 simulation)
   * - track ID: Geant4 track ID of the ionizing particle
   * - energy: amount of energy released by ionization (from Geant4 simulation)
   * - electrons: amount of electrons reaching the readout channel
   *
   * Note the different definition of the electrons respect to the rest:
   * it describes the electrons at the anode _after_ the drifting occurred,
   * while all the other quantities can be related to the moment the ionization
   * happened.
   *
   * The number of electrons typically includes inefficiencies and physics
   * effects that reduce and spread the electrons. In the simulation,
   * this yields a fractional number of electrons.
   *
   * Each IDE is also typically associated with a time (TDC) count, that is
   * the time at which the ionized electrons reached the readout channel, in
   * electronic ticks, as opposed as the time when ionization occurred.
   * The latter is not stored.
   *
   * At the time of writing this documentation (LArSoft 6.4.0), IDEs are
   * computed in larg4::LArVoxelReadout.
   * The energy and track ID come directly from Geant4 simulation.
   * The position is the mid point of the Geant4 step that produced ionization.
   * The electrons are
   *
   * 1. converted from that same energy (using a fundamental conversion factor
   *    stored in `larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h`)
   * 2. applied recombination effect by larg4::IonizationAndScintillation::Reset()
   * 3. applied attenuation and diffusion in
   *    larg4::LArVoxelReadout::DriftIonizationElectrons()
   *
   * The latter also assembles the sim::IDE objects to be stored into
   * sim::SimChannel.
   *
   */
  struct IDE {

    /// Type of track ID (the value comes from Geant4)
    typedef int TrackID_t;

    /// Default constructor (sets "bogus" values)
    IDE();

    /// Constructor: copies an IDE, and applies the specified offset to track ID
    IDE(IDE const& ide, int offset);

    /// Constructor: sets all data members
    IDE(TrackID_t tid,
        float nel,
        float e,
        float xpos,
        float ypos,
        float zpos,
        TrackID_t gid = util::kBogusI)
      : trackID(tid), numElectrons(nel), energy(e), x(xpos), y(ypos), z(zpos), origTrackID(gid)
    {}

    TrackID_t trackID;  ///< Geant4 supplied track ID
    float numElectrons; ///< number of electrons at the readout for this track ID and time
    float energy;       ///< energy deposited by ionization by this track ID and time [MeV]
    float x;            ///< x position of ionization [cm]
    float y;            ///< y position of ionization [cm]
    float z;            ///< z position of ionization [cm]
    TrackID_t
      origTrackID; ///< Geant4 supplied track ID (remains true trackID even for shower secondaries/tertiaries etc)
  };               // struct IDE

  /// List of energy deposits at the same time (on this channel)
  typedef std::pair<unsigned short, std::vector<sim::IDE>> TDCIDE;

  /**
   * @brief Energy deposited on a readout channel by simulated tracks
   *
   * This class stores the list of all energies deposited on a readout channel.
   * The number of electrons is stored as well.
   *
   * The information is organized by time: it is divided by TDC ticks, and
   * each TDC tick where some energy was deposited appears in a separate entry,
   * while the quiet TDC ticks are omitted.
   * For each TDC, the information is stored as a list of energy deposits;
   * each deposit comes from a single Geant4 track and stores the location where
   * the ionization happened according to the simulation (see `sim::IDE` class).
   *
   * Note that there can be multiple energy deposit records (that is `sim::IDE`)
   * for a single track in a single TDC tick.
   */
  class SimChannel {
  public:
    /// Type for TDC tick used in the internal representation
    typedef TDCIDE::first_type StoredTDC_t;

    /// Type of list of energy deposits for each TDC with signal
    typedef std::vector<TDCIDE> TDCIDEs_t;

  private:
    /// Temporary opaque replacement for `std::span<TDCIDE>`.
    /// @todo Directly replace with `std::span<TDCIDE>` with C++20.
    struct TDCIDEspan_t {
      typedef TDCIDEs_t::const_iterator iterator;
      iterator b, e;
      iterator begin() const noexcept { return b; }
      iterator end() const noexcept { return e; }
    };

    raw::ChannelID_t fChannel; ///< readout channel where electrons are collected
    TDCIDEs_t fTDCIDEs;        ///< list of energy deposits for each TDC with signal

  public:
    // Default constructor
    SimChannel();

    /// Type for TDC tick used in the interface
    /// (different type than raw::TDCtick_t! and from internal representation!
    ///  but same meaning!)
    typedef unsigned int TDC_t;

    /// Type of track ID (the value comes from Geant4)
    typedef IDE::TrackID_t TrackID_t;

    /// Type of iterator to a `sim::TDCIDE`.
    using TDCIDEIter_t = TDCIDEs_t::const_iterator;

    /// Constructor: immediately sets the channel number
    explicit SimChannel(raw::ChannelID_t channel);

    /**
     * @brief Add ionization electrons and energy to this channel
     * @param trackID ID of simulated track depositing this energy (from Geant4)
     * @param tdc TDC tick when this deposit was collected
     * @param numberElectrons electrons created at this point by this track
     * @param xyz coordinates of original location of ionization (3D array) [cm]
     * @param energy energy deposited at this point by this track [MeV]
     *
     * The number of electrons can be fractional because of simulated
     * efficiency and physics effects.
     */
    void AddIonizationElectrons(TrackID_t trackID,
                                TDC_t tdc,
                                double numberElectrons,
                                double const* xyz,
                                double energy,
                                TrackID_t origTrackID = util::kBogusI);

    /// Returns the readout channel this object describes
    raw::ChannelID_t Channel() const;

    /**
     * @brief Return all the recorded energy deposition within a time interval
     * @param startTDC TDC tick opening the time window
     * @param endTDC TDC tick closing the time window (included in the interval)
     * @return a collection of energy deposit information from all tracks
     *
     * This method returns the energy deposited on this channel by each track
     * ID active in the specified TDC time interval.
     *
     * Each entry pertains a single track ID. For each entry, all energy
     * deposit information is merged into a single record. It includes:
     * * energy and number of electrons, as the integral in the time interval
     * * position, as average weighted by the number of electrons
     * * the ID of the track depositing this energy
     *
     * Entries are sorted by track ID number.
     */
    std::vector<sim::IDE> TrackIDsAndEnergies(TDC_t startTDC, TDC_t endTDC) const;

    /**
     * @brief Returns all the deposited energy information as stored
     * @return all the deposited energy information as stored in the object
     *
     * The returned list is organized in pairs. Each pair contains all
     * ionization information in a single TDC tick (collection of `sim::IDE`),
     * and the number of that tick. The information is sorted by increasing TDC
     * tick.
     *
     * See the class description for the details of the ionization information
     * content.
     */
    TDCIDEs_t const& TDCIDEMap() const;

    /**
     * @brief Returns iterators to first-and past-last `TDCIDE` in tick interval.
     * @param from tick of the time interval start
     * @param to tick of the time interval end
     * @return iterator span between first `TDCIDE` no earlier than `from` and than `to`
     *
     * The returned value supports data members:
     *  * `begin()`: iterator to the first `sim::TDCIDE` at or after tick `from`
     *  * `end()`: iterator to the first `sim::TDCIDE` at or after tick `to`
     *
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * double chargeBetween(sim::SimChannel const& sc, int start, int ticks) {
     *   double charge = 0.0;
     *   for (sim::TDCIDE const& TDCIDE: sc.IDEsBetween(start, start + ticks)) {
     *     for (sim::IDE const& ide: TDCIDE.second)
     *       charge += ide.numElectrons;
     *   }
     *   return charge;
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * All `sim::TDCIDE` in the returned span have TPC tick at `from` or after,
     * and strictly before `to`.
     */
    auto IDEsBetween(TDC_t from, TDC_t to) const;

    /**
     * @brief Returns the first `sim::TDCIDE` at or after the specified `tick`.
     * @param tick a TDC tick number
     * @return pointer to the first `sim::TDCIDE`, or `nullptr` if none
     *
     * The returned `sim::TDCIDE` is the first at TDC value `tick` or after it.
     * It is therefore not guaranteed that the returned IDE are _exactly_ at
     * `tick`.
     * If `tick` is larger than the last tick in the object, `nullptr` is
     * returned.
     */
    sim::TDCIDE const* findIDEatTime(TDC_t tick) const;

    /// Returns the total number of ionization electrons on this channel in the specified TDC
    double Charge(TDC_t tdc) const;

    /// Returns the total energy on this channel in the specified TDC [MeV]
    double Energy(TDC_t tdc) const;

    /**
     * @brief Returns energies collected for each track within a time interval
     * @param startTDC TDC tick opening the time window
     * @param endTDC TDC tick closing the time window (included in the interval)
     * @return a collection of energy and fraction from each track in interval
     * @see TrackIDsAndEnergies()
     *
     * This method returns the energy deposited on this channel by each track
     * ID active in the specified TDC time interval.
     *
     * Each entry pertains a single track ID. For each entry, all energy
     * deposit information is merged into a single record. It includes:
     * * energy of the track, as the integral in the time interval [MeV]
     * * energy fraction respect to the total (see below)
     * * the ID of the track depositing this energy
     *
     * The energy fraction is the energy deposited by the track on this channel
     * in the specified time interval, divided by the total of the energy
     * deposited by all tracks on this channel in that same time interval.
     *
     * Entries are sorted by track ID number.
     */
    std::vector<sim::TrackIDE> TrackIDEs(TDC_t startTDC, TDC_t endTDC) const;

    /// Comparison: sorts by channel ID
    bool operator<(const SimChannel& other) const;

    /// Comparison: true if SimChannels have the same channel ID
    bool operator==(const SimChannel& other) const;

    /**
     * @brief Merges the deposits from another channel into this one
     * @param channel the sim::SimChannel holding information to be merged
     * @param offset track ID offset for the merge
     * @return range of the IDs of the added tracks
     *
     * The information from the specified simulated channel is added to the
     * current one.
     * This is achieved by appending the energy deposit information (`sim::IDE`)
     * at each TDC tick from the merged channel to the list of existing energy
     * deposits for that TDC tick.
     *
     * In addition, the track IDs of the merged channel are added an offset,
     * so that they can be distinguished from the existing ones.
     * This is useful when simulating tracks with multiple Geant4 runs. Geant4
     * will reuse track IDs on each run, and using the highest number of track
     * ID from one run as the offset for the next avoids track ID collisions.
     * Note however that this function does not perform any collision check, and
     * it is caller's duty to ensure that the offset is properly large.
     * The return value is a pair including the lowest and the largest track IDs
     * added to this channel, equivalent to the lowest and the highest track IDs
     * present in the merged channel, both augmented by the applied offset.
     *
     * The channel number of the merged channel is ignored.
     */
    std::pair<TrackID_t, TrackID_t> MergeSimChannel(const SimChannel& channel, int offset);

    /**
     * @brief Dumps the full content of the SimChannel into a stream
     * @tparam Stream an ostream-like stream object
     * @param out the stream to send the information into
     * @param indent indentation of the lines (default: none)
     * @param first_indent indentation for the first line (default: as `indent`)
     */
    template <typename Stream>
    void Dump(Stream&& out, std::string indent, std::string first_indent) const;

    /// Documentation at `Dump(Stream&&, std::string, std::string) const`.
    template <typename Stream>
    void Dump(Stream&& out, std::string indent = "") const
    {
      Dump(std::forward<Stream>(out), indent, indent);
    }

  private:
    /// Comparison functor, sorts by increasing TDCtick value
    struct CompareByTDC;

    /// Return the iterator to the first TDCIDE not earlier than tdc
    TDCIDEs_t::iterator findClosestTDCIDE(StoredTDC_t tdc);

    /// Return the (constant) iterator to the first TDCIDE not earlier than tdc
    TDCIDEs_t::const_iterator findClosestTDCIDE(StoredTDC_t tdc) const;
    /// @}
  };

} // namespace sim

inline bool sim::SimChannel::operator<(const sim::SimChannel& other) const
{
  return fChannel < other.Channel();
}
inline bool sim::SimChannel::operator==(const sim::SimChannel& other) const
{
  return fChannel == other.Channel();
}
inline sim::SimChannel::TDCIDEs_t const& sim::SimChannel::TDCIDEMap() const
{
  return fTDCIDEs;
}
inline auto sim::SimChannel::IDEsBetween(TDC_t from, TDC_t to) const
{
  return TDCIDEspan_t{findClosestTDCIDE(from), findClosestTDCIDE(to)};
}
inline sim::TDCIDE const* sim::SimChannel::findIDEatTime(TDC_t tick) const
{
  auto const it = findClosestTDCIDE(tick);
  return (it == TDCIDEMap().end()) ? nullptr : &*it;
}
inline raw::ChannelID_t sim::SimChannel::Channel() const
{
  return fChannel;
}

// -----------------------------------------------------------------------------
// ---  template implementation
// ---
template <class Stream>
void sim::SimChannel::Dump(Stream&& out, std::string indent, std::string first_indent) const
{
  out << first_indent << "channel #" << Channel() << " read " << fTDCIDEs.size() << " TDCs:\n";
  double channel_energy = 0., channel_charge = 0.;
  for (const auto& TDCinfo : fTDCIDEs) {
    auto const tdc = TDCinfo.first;
    out << indent << "  TDC #" << tdc << " with " << TDCinfo.second.size() << " IDEs\n";
    double tdc_energy = 0., tdc_charge = 0.;
    for (const sim::IDE& ide : TDCinfo.second) {
      out << indent << "    (" << ide.x << ", " << ide.y << ", " << ide.z << ") "
          << ide.numElectrons << " electrons, " << ide.energy << " MeV (trkID=" << ide.trackID
          << ")\n";
      tdc_energy += ide.energy;
      tdc_charge += ide.numElectrons;
    } // for IDEs
    out << indent << "    => TDC #" << tdc << " CH #" << Channel() << " collected " << tdc_charge
        << " electrons and " << tdc_energy << " MeV\n";
    channel_energy += tdc_energy;
    channel_charge += tdc_charge;
  } // for TDCs
  out << indent << "  => channel #" << Channel() << " collected " << channel_charge
      << " electrons and " << channel_energy << " MeV\n";
} // sim::SimChannel::Dump<>()

#endif // LARDATAOBJ_SIMULATION_SIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
