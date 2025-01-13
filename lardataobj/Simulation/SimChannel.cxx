///
/// \file  Simulation/SimChannel.cxx
///
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#include <algorithm> // std::lower_bound(), std::max()
#include <limits>    // std::numeric_limits
#include <map>
#include <stdexcept>
#include <utility>

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/sim.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace sim {

  //-------------------------------------------------
  IDE::IDE()
    : trackID(util::kBogusI)
    , numElectrons(util::kBogusD)
    , energy(util::kBogusD)
    , x(util::kBogusD)
    , y(util::kBogusD)
    , z(util::kBogusD)
    , origTrackID(util::kBogusI)
  {}

  //-------------------------------------------------
  IDE::IDE(sim::IDE const& ide, int offset) : IDE(ide)
  {

    trackID = trackID >= 0 ? trackID + offset : trackID - offset;
    origTrackID = origTrackID >= 0 ? origTrackID + offset : origTrackID - offset;
  }

  // Default constructor
  //-------------------------------------------------
  SimChannel::SimChannel() : fChannel(raw::InvalidChannelID) {}

  //-------------------------------------------------
  SimChannel::SimChannel(raw::ChannelID_t channel) : fChannel(channel) {}

  //-------------------------------------------------
  void SimChannel::AddIonizationElectrons(TrackID_t trackID,
                                          TDC_t tdc,
                                          double numberElectrons,
                                          double const* xyz,
                                          double energy,
                                          TrackID_t origTrackID)
  {
    // look at the collection to see if the current TDC already
    // exists, if not, add it, if so, just add a new track id to the
    // vector, or update the information if track is already present

    // no electrons? no energy? no good!
    if ((numberElectrons < std::numeric_limits<double>::epsilon()) ||
        (energy <= std::numeric_limits<double>::epsilon())) {
      // will throw
      MF_LOG_ERROR("SimChannel") << "AddIonizationElectrons() trying to add to TDC #" << tdc << " "
                                 << numberElectrons << " electrons with " << energy
                                 << " MeV of energy from track ID=" << trackID;
      return;
    } // if no energy or no electrons

    auto itr = findClosestTDCIDE(tdc);

    // check if this tdc value is in the vector, it is possible that
    // the lower bound is different from the given tdc, in which case
    // we need to add something for that tdc
    if (itr == fTDCIDEs.end() || itr->first != tdc) {
      std::vector<sim::IDE> idelist;
      idelist.emplace_back(trackID, numberElectrons, energy, xyz[0], xyz[1], xyz[2], origTrackID);
      fTDCIDEs.emplace(itr, tdc, std::move(idelist));
    }
    else { // we have that TDC already; itr points to it

      // loop over the IDE vector for this tdc and add the electrons
      // to the entry with the same G4 track id
      for (auto& ide : itr->second) {

        if (ide.origTrackID != origTrackID) continue;

        // make a weighted average for the location information
        double weight = ide.numElectrons + numberElectrons;
        ide.x = (ide.x * ide.numElectrons + xyz[0] * numberElectrons) / weight;
        ide.y = (ide.y * ide.numElectrons + xyz[1] * numberElectrons) / weight;
        ide.z = (ide.z * ide.numElectrons + xyz[2] * numberElectrons) / weight;
        ide.numElectrons = weight;
        ide.energy = ide.energy + energy;

        // found the track id we wanted, so return;
        return;
      } // for

      // if we never found the track id, then this is the first instance of
      // the track id for this tdc, so add ide to the vector
      itr->second.emplace_back(
        trackID, numberElectrons, energy, xyz[0], xyz[1], xyz[2], origTrackID);

    } // if new TDC ... else

  } // SimChannel::AddIonizationElectrons()

  //-------------------------------------------------
  double SimChannel::Charge(TDC_t tdc) const
  {
    double charge = 0.;

    auto itr = findClosestTDCIDE(tdc);

    // check to see if this tdc value is in the map
    if (itr != fTDCIDEs.end() && itr->first == tdc) {

      // loop over the list for this tdc value and add up
      // the total number of electrons
      for (auto ide : itr->second) {
        charge += ide.numElectrons;
      } // end loop over sim::IDE for this tdc

    } // end if this tdc is represented in the map

    return charge;
  }

  //-------------------------------------------------
  double SimChannel::Energy(TDC_t tdc) const
  {
    double energy = 0.;

    auto itr = findClosestTDCIDE(tdc);

    // check to see if this tdc value is in the map
    if (itr != fTDCIDEs.end() && itr->first == tdc) {

      // loop over the list for this tdc value and add up
      // the total number of electrons
      for (auto ide : itr->second) {
        energy += ide.energy;
      } // end loop over sim::IDE for this tdc

    } // end if this tdc is represented in the map

    return energy;
  }

  //-----------------------------------------------------------------------
  // the start and end tdc values are assumed to be inclusive
  std::vector<sim::IDE> SimChannel::TrackIDsAndEnergies(TDC_t startTDC, TDC_t endTDC) const
  {
    // make a map of track ID values to sim::IDE objects

    if (startTDC > endTDC) {
      mf::LogWarning("SimChannel")
        << "requested tdc range is bogus: " << startTDC << " " << endTDC << " return empty vector";
      return {}; // returns an empty vector
    }

    std::map<TrackID_t, sim::IDE> idToIDE;

    //find the lower bound for this tdc and then iterate from there
    auto itr = findClosestTDCIDE(startTDC);

    while (itr != fTDCIDEs.end()) {

      // check the tdc value for the iterator, break the loop if we
      // are outside the range
      if (itr->first > endTDC) break;

      // grab the vector of IDEs for this tdc
      auto const& idelist = itr->second;
      // now loop over them and add their content to the map
      for (auto const& ide : idelist) {
        auto itTrkIDE = idToIDE.find(ide.trackID);
        if (itTrkIDE != idToIDE.end()) {
          // the IDE we are going to update:
          sim::IDE& trackIDE = itTrkIDE->second;

          double const nel1 = trackIDE.numElectrons;
          double const nel2 = ide.numElectrons;
          double const en1 = trackIDE.energy;
          double const en2 = ide.energy;
          double const energy = en1 + en2;
          double const weight = nel1 + nel2;

          // make a weighted average for the location information
          trackIDE.x = (ide.x * nel2 + trackIDE.x * nel1) / weight;
          trackIDE.y = (ide.y * nel2 + trackIDE.y * nel1) / weight;
          trackIDE.z = (ide.z * nel2 + trackIDE.z * nel1) / weight;
          trackIDE.numElectrons = weight;
          trackIDE.energy = energy;
        } // end if the track id for this one is found
        else {
          idToIDE[ide.trackID] = sim::IDE(ide);
        }
      } // end loop over vector

      ++itr;
    } // end loop over tdc values

    // now fill the vector with the ides from the map
    std::vector<sim::IDE> ides;
    ides.reserve(idToIDE.size());
    for (auto const& itr : idToIDE) {
      ides.push_back(itr.second);
    }

    return ides;
  }

  //-----------------------------------------------------------------------
  // the start and end tdc values are assumed to be inclusive
  std::vector<sim::TrackIDE> SimChannel::TrackIDEs(TDC_t startTDC, TDC_t endTDC) const
  {

    std::vector<sim::TrackIDE> trackIDEs;

    if (startTDC > endTDC) {
      mf::LogWarning("SimChannel::TrackIDEs")
        << "requested tdc range is bogus: " << startTDC << " " << endTDC << " return empty vector";
      return trackIDEs;
    }

    double totalE = 0.;
    std::vector<sim::IDE> const ides = TrackIDsAndEnergies(startTDC, endTDC);
    for (auto const& ide : ides)
      totalE += ide.energy;

    // protect against a divide by zero below
    if (totalE < 1.e-5) totalE = 1.;

    // loop over the entries in the map and fill the input vectors
    for (auto const& ide : ides) {
      if (ide.trackID == sim::NoParticleId) continue;
      trackIDEs.emplace_back(
        ide.trackID, ide.energy / totalE, ide.energy, ide.numElectrons, ide.origTrackID);
    }

    return trackIDEs;
  }

  //-----------------------------------------------------------------------
  // Merge the collection of IDEs from one sim channel to another.
  // Requires an agreed upon offset for G4 trackID
  std::pair<SimChannel::TrackID_t, SimChannel::TrackID_t> SimChannel::MergeSimChannel(
    SimChannel const& channel,
    int offset)
  {
    if (this->Channel() != channel.Channel())
      throw std::runtime_error("ERROR SimChannel Merge: Trying to merge different channels!");

    std::pair<TrackID_t, TrackID_t> range_trackID(std::numeric_limits<int>::max(),
                                                  std::numeric_limits<int>::min());

    for (auto const& itr : channel.TDCIDEMap()) {

      auto tdc = itr.first;
      auto const& ides = itr.second;

      // find the entry from this SimChannel corresponding to the tdc from the other
      auto itrthis = findClosestTDCIDE(tdc);

      // if not found, insert at the location where it should be
      if (itrthis == fTDCIDEs.end() || itrthis->first != tdc) {
        itrthis = fTDCIDEs.insert(itrthis, {tdc, std::vector<sim::IDE>()});
      }

      // pick the IDE list we have to fill
      std::vector<sim::IDE>* curIDEVec = &(itrthis->second);

      for (auto const& ide : ides) {
        curIDEVec->emplace_back(ide, offset);
        auto tid = std::abs(ide.trackID) + offset;
        if (tid < range_trackID.first) range_trackID.first = tid;
        if (tid > range_trackID.second) range_trackID.second = tid;
      } //end loop over IDEs

    } //end loop over TDCIDEMap

    return range_trackID;
  }

  //-------------------------------------------------
  struct SimChannel::CompareByTDC {

    bool operator()(TDCIDE const& a, TDCIDE const& b) const { return a.first < b.first; }

    bool operator()(StoredTDC_t a_tdc, TDCIDE const& b) const { return a_tdc < b.first; }

    bool operator()(TDCIDE const& a, StoredTDC_t b_tdc) const { return a.first < b_tdc; }

  }; // struct CompareByTDC

  SimChannel::TDCIDEs_t::iterator SimChannel::findClosestTDCIDE(StoredTDC_t tdc)
  {
    return std::lower_bound(fTDCIDEs.begin(), fTDCIDEs.end(), tdc, CompareByTDC());
  }

  SimChannel::TDCIDEs_t::const_iterator SimChannel::findClosestTDCIDE(StoredTDC_t tdc) const
  {
    return std::lower_bound(fTDCIDEs.begin(), fTDCIDEs.end(), tdc, CompareByTDC());
  }

  //-------------------------------------------------

}
