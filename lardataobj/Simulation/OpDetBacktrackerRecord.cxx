///
/// \file  Simulation/OpDetBacktrackerRecord.cxx
///
///
/// \author  jason.stock@mines.sdsmt.edu
//
//   based on the SimChannel object by seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#include <algorithm> // std::lower_bound(), std::max()
#include <limits>    // std::numeric_limits
#include <map>
#include <stdexcept>
#include <utility>

#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/Simulation/OpDetBacktrackerRecord.h"
#include "lardataobj/Simulation/sim.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
namespace sim {

  //-------------------------------------------------
  SDP::SDP() //SDP stands for Scintillation Deposited Photons
    : trackID(util::kBogusI)
    , numPhotons(util::kBogusD)
    , energy(util::kBogusD)
    , x(util::kBogusD)
    , y(util::kBogusD)
    , z(util::kBogusD)
  {}

  //-------------------------------------------------
  SDP::SDP(sim::SDP const& sdp, int offset) : SDP(sdp)
  {
    trackID += offset;
  }

  // Default constructor
  //-------------------------------------------------
  OpDetBacktrackerRecord::OpDetBacktrackerRecord()
    : iOpDetNum(
        -1) //set an impossible channel number in the case where this is called without an opticalchannel.
  //The reason for doing this is to follow the structure of SimChannel, which uses kBogusChannel
  {}

  //-------------------------------------------------
  OpDetBacktrackerRecord::OpDetBacktrackerRecord(int detNum) : iOpDetNum(detNum) {}

  //-------------------------------------------------
  void OpDetBacktrackerRecord::AddScintillationPhotons(TrackID_t trackID,
                                                       timePDclock_t iTimePDclock,
                                                       double numberPhotons,
                                                       double const* xyz,
                                                       double energy)
  {
    /**
     * The iTimePDclock from OpFastScintillation (where OpDetBacktrackerRecords originate) is done with CLHEP::ns (units of nanoseconds).
     */
    // look at the collection to see if the current iTimePDclock already
    // exists, if not, add it, if so, just add a new track id to the
    // vector, or update the information if track is already present

    // no photons? no good!
    if ((numberPhotons < std::numeric_limits<double>::epsilon()) ||
        (energy <= std::numeric_limits<double>::epsilon())) {
      // will throw
      MF_LOG_ERROR("OpDetBacktrackerRecord")
        << "AddTrackPhotons() trying to add to iTimePDclock #" << iTimePDclock << " "
        << numberPhotons << " photons with " << energy
        << " MeV of energy from track ID=" << trackID;
      return;
    } // if no photons

    int rounded_time = std::round(iTimePDclock);

    // auto itr = findClosestTimePDclockSDP(rounded_time);
    auto itr = timePDclockSDPs.lower_bound(rounded_time);

    // check if this iTimePDclock value is in the vector, it is possible that
    // the lower bound is different from the given TimePDclock, in which case
    // we need to add something for that TimePDclock
    if (itr == timePDclockSDPs.end() || timePDclockSDPs.key_comp()(rounded_time, itr->first)) {
      timePDclockSDPs.insert(
        {rounded_time,
         std::vector<sim::SDP>{sim::SDP(trackID, numberPhotons, energy, xyz[0], xyz[1], xyz[2])}});
    }
    else { // we have that iTimePDclock already; itr points to it

      // loop over the SDP vector for this iTimePDclock and add the electrons
      // to the entry with the same track id
      for (auto& sdp : itr->second) {

        if (sdp.trackID != trackID) continue;

        // make a weighted average for the location information
        double weight = sdp.numPhotons + numberPhotons;
        sdp.x = (sdp.x * sdp.numPhotons + xyz[0] * numberPhotons) / weight;
        sdp.y = (sdp.y * sdp.numPhotons + xyz[1] * numberPhotons) / weight;
        sdp.z = (sdp.z * sdp.numPhotons + xyz[2] * numberPhotons) / weight;
        sdp.numPhotons = weight;
        sdp.energy = sdp.energy + energy;

        // found the track id we wanted, so return;
        return;
      } // for

      // if we never found the track id, then this is the first instance of
      // the track id for this TimePDclock, so add sdp to the vector
      itr->second.emplace_back(trackID, numberPhotons, energy, xyz[0], xyz[1], xyz[2]);

    } // end of else "We have that iTimePDclock already"

  } // OpDetBacktrackerRecord::AddIonizationElectrons()

  //-------------------------------------------------
  double OpDetBacktrackerRecord::Photons(timePDclock_t iTimePDclock) const //Number of photons
  {
    double numPhotons = 0.;

    // auto itr = findClosestTimePDclockSDP(iTimePDclock);
    auto itr = timePDclockSDPs.lower_bound(iTimePDclock);
    // check to see if this iTimePDclock value is in the map
    if (itr != timePDclockSDPs.end() && itr->first == iTimePDclock) {

      // loop over the list for this iTimePDclock value and add up
      // the total number of electrons
      for (auto sdp : itr->second) {
        numPhotons += sdp.numPhotons;
      } // end loop over sim::SDP for this TimePDclock

    } // end if this iTimePDclock is represented in the map

    return numPhotons;
  }

  //-------------------------------------------------
  //Energy in each object is the energy deposited by the track step that left the photons.
  double OpDetBacktrackerRecord::Energy(timePDclock_t iTimePDclock) const
  {
    double energy = 0.;

    // auto itr = findClosestTimePDclockSDP(iTimePDclock);
    auto itr = timePDclockSDPs.lower_bound(iTimePDclock);

    // check to see if this iTimePDclock value is in the map
    if (itr != timePDclockSDPs.end() && itr->first == iTimePDclock) {

      // loop over the list for this iTimePDclock value and add up
      // the total number of photons
      for (auto sdp : itr->second) {
        energy += sdp.energy;
      } // end loop over sim::SDP for this TimePDclock

    } // end if this iTimePDclock is represented in the map

    return energy;
  }

  //-----------------------------------------------------------------------
  // the start and end iTimePDclock values are assumed to be inclusive
  std::vector<sim::SDP> OpDetBacktrackerRecord::TrackIDsAndEnergies(
    timePDclock_t startTimePDclock,
    timePDclock_t endTimePDclock) const
  {
    // make a map of track ID values to sim::SDP objects

    if (startTimePDclock > endTimePDclock) {
      mf::LogWarning("OpDetBacktrackerRecord")
        << "requested TimePDclock range is bogus: " << startTimePDclock << " " << endTimePDclock
        << " return empty vector";
      return {}; // returns an empty vector
    }

    std::map<TrackID_t, sim::SDP> idToSDP;

    //find the lower bound for this iTimePDclock and then iterate from there
    // auto itr = findClosestTimePDclockSDP(startTimePDclock);
    auto itr = timePDclockSDPs.lower_bound(startTimePDclock);

    while (itr != timePDclockSDPs.end()) {

      // check the iTimePDclock value for the iterator, break the loop if we
      // are outside the range
      if (itr->first > endTimePDclock) break;

      // grab the vector of SDPs for this TimePDclock
      auto const& sdplist = itr->second;
      // now loop over them and add their content to the map
      for (auto const& sdp : sdplist) {
        auto itTrkSDP = idToSDP.find(sdp.trackID);
        if (itTrkSDP != idToSDP.end()) {
          // the SDP we are going to update:
          sim::SDP& trackSDP = itTrkSDP->second;

          double const nPh1 = trackSDP.numPhotons;
          double const nPh2 = sdp.numPhotons;
          double const weight = nPh1 + nPh2;

          // make a weighted average for the location information
          trackSDP.x = (sdp.x * nPh2 + trackSDP.x * nPh1) / weight;
          trackSDP.y = (sdp.y * nPh2 + trackSDP.y * nPh1) / weight;
          trackSDP.z = (sdp.z * nPh2 + trackSDP.z * nPh1) / weight;
          trackSDP.numPhotons = weight;
        } // end if the track id for this one is found
        else {
          idToSDP[sdp.trackID] = sim::SDP(sdp);
        }
      } // end loop over vector

      ++itr;
    } // end loop over iTimePDclock values

    // now fill the vector with the sdps from the map
    std::vector<sim::SDP> sdps;
    sdps.reserve(idToSDP.size());
    for (auto const& itr : idToSDP) {
      sdps.push_back(itr.second);
    }

    return sdps;
  }

  //-----------------------------------------------------------------------
  // the start and end iTimePDclock values are assumed to be inclusive
  std::vector<sim::TrackSDP> OpDetBacktrackerRecord::TrackSDPs(timePDclock_t startTimePDclock,
                                                               timePDclock_t endTimePDclock) const
  {

    std::vector<sim::TrackSDP> trackSDPs;

    if (startTimePDclock > endTimePDclock) {
      mf::LogWarning("OpDetBacktrackerRecord::TrackSDPs")
        << "requested iTimePDclock range is bogus: " << startTimePDclock << " " << endTimePDclock
        << " return empty vector";
      return trackSDPs;
    }

    double totalPhotons = 0.;
    std::vector<sim::SDP> const sdps = TrackIDsAndEnergies(startTimePDclock, endTimePDclock);
    for (auto const& sdp : sdps)
      totalPhotons += sdp.numPhotons;

    // protect against a divide by zero below
    if (totalPhotons < 1.e-5) totalPhotons = 1.;

    // loop over the entries in the map and fill the input vectors
    for (auto const& sdp : sdps) {
      if (sdp.trackID == sim::NoParticleId) continue;
      trackSDPs.emplace_back(sdp.trackID, sdp.numPhotons / totalPhotons, sdp.numPhotons);
    }

    return trackSDPs;
  }

  //-----------------------------------------------------------------------
  // Merge the collection of SDPs from one sim channel to another.
  // Requires an agreed upon offset for G4 trackID
  std::pair<OpDetBacktrackerRecord::TrackID_t, OpDetBacktrackerRecord::TrackID_t>
  OpDetBacktrackerRecord::MergeOpDetBacktrackerRecord(OpDetBacktrackerRecord const& channel,
                                                      int offset)
  {
    if (this->OpDetNum() != channel.OpDetNum())
      throw std::runtime_error(
        "ERROR OpDetBacktrackerRecord Merge: Trying to merge different channels!");

    std::pair<TrackID_t, TrackID_t> range_trackID(std::numeric_limits<int>::max(),
                                                  std::numeric_limits<int>::min());

    for (auto const& [iTimePDclock, sdps] : channel.timePDclockSDPsMap()) {

      // find the entry from this OpDetBacktrackerRecord corresponding to the iTimePDclock from the other
      auto itrthis = timePDclockSDPs.lower_bound(iTimePDclock);
      // pick which SDP list we have to fill: new one or existing one
      std::vector<sim::SDP>* curSDPVec;
      if (itrthis == timePDclockSDPs.end() || itrthis->first != iTimePDclock) {
        curSDPVec =
          &(timePDclockSDPs.insert({iTimePDclock, std::vector<sim::SDP>()}).first->second);
      }
      else
        curSDPVec = &(itrthis->second);

      for (auto const& sdp : sdps) {
        curSDPVec->emplace_back(sdp, offset);
        if (sdp.trackID + offset < range_trackID.first) range_trackID.first = sdp.trackID + offset;
        if (sdp.trackID + offset > range_trackID.second)
          range_trackID.second = sdp.trackID + offset;
      } //end loop over SDPs

    } //end loop over TimePDclockSDPMap

    return range_trackID;
  }
}
