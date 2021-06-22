/**
 * @file    Wire_test.cc
 * @brief   Simple test on a recob::Wire object
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20150113
 * @version 1.0
 *
 * This test simply creates recob::Wire objects and verifies that the values it
 * can access are the right ones.
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * version 1.0: ~1.5" (debug mode)
 */

// C/C++ standard library
#include <algorithm> // std::equal()


// Boost libraries
/*
 * Boost Magic: define the name of the module;
 * and do that before the inclusion of Boost unit test headers
 * because it will change what they provide.
 * Among the those, there is a main() function and some wrapping catching
 * unhandled exceptions and considering them test failures, and probably more.
 * This also makes fairly complicate to receive parameters from the command line
 * (for example, a random seed).
 */
#define BOOST_TEST_MODULE ( wire_test )
#include "cetlib/quiet_unit_test.hpp" // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK()

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t
#include "lardataobj/Utilities/sparse_vector.h"
#include "lardataobj/RecoBase/Wire.h"



//------------------------------------------------------------------------------
//--- Test code
//


void CheckWire(
  recob::Wire const& wire,
  recob::Wire::RegionsOfInterest_t const& sigROIlist,
  raw::ChannelID_t channel,
  geo::View_t view
) {

  // verify that the values are as expected
  // - channel ID
  BOOST_CHECK_EQUAL(wire.Channel(), channel);

  // - view
  BOOST_CHECK_EQUAL(wire.View(), view);

  // - region of interest
  BOOST_CHECK_EQUAL(wire.NSignal(), sigROIlist.size());

  recob::Wire::RegionsOfInterest_t const& wireROI = wire.SignalROI();
  BOOST_CHECK_EQUAL(wireROI.n_ranges(), sigROIlist.n_ranges());

  unsigned int index = 0;
  for (auto sample: wireROI) {
    BOOST_CHECK_EQUAL(sample, sigROIlist[index++]);
  }

  // - other elements of interface
  auto const& wire_signal = wire.Signal();
  BOOST_CHECK
    (std::equal(wire_signal.begin(), wire_signal.end(), sigROIlist.cbegin()));

} // CheckWire()


void WireTestDefaultConstructor() {

  //
  // Part I: initialization of wire inputs
  //
  // these are the values expected for a default-constructed wire
  raw::ChannelID_t channel = raw::InvalidChannelID;
  geo::View_t view = geo::kUnknown;
  recob::Wire::RegionsOfInterest_t sigROIlist;

  //
  // Part II: default constructor
  //
  // step II.1: create a wire with the default constructor
  recob::Wire wire;


  // step II.2: verify that the values are as expected
  CheckWire(wire, sigROIlist, channel, view);

} // WireTestDefaultConstructor()


void WireTestCustomConstructors() {

  //
  // Part I: initialization of wire inputs
  //
  raw::ChannelID_t channel = 12;
  geo::View_t view = geo::kV;

  recob::Wire::RegionsOfInterest_t sigROIlist(20);
  sigROIlist.add_range
    (5, recob::Wire::RegionsOfInterest_t::vector_t({ 5., 6., 7. }));
  sigROIlist.add_range
    (11, recob::Wire::RegionsOfInterest_t::vector_t({ 11., 12., 13., 14. }));

  // this is not a recob::Wire test, but I want to make sure anyway...
  BOOST_CHECK_EQUAL(sigROIlist.size(), 20U);
  BOOST_CHECK_EQUAL(sigROIlist.n_ranges(), 2U);
  size_t index = 0;
  for (auto sample: sigROIlist) {
    if(sample != 0) {BOOST_CHECK((sample == (float) index) );}
    ++index;
  } // for


  //
  // Part II: constructor with signal copy
  //
  // step II.1: create a wire with the signal-copying constructor
  recob::Wire wire1(sigROIlist, channel, view);


  // step II.2: verify that the values are as expected
  CheckWire(wire1, sigROIlist, channel, view);


  //
  // Part III: constructor with signal move
  //
  // step III.1: create a wire with the signal-copying constructor
  recob::Wire::RegionsOfInterest_t sigROIlistCopy(sigROIlist); // need a copy for check
  recob::Wire wire2(std::move(sigROIlistCopy), channel, view);

  // step III.2: verify that the values are as expected
  CheckWire(wire2, sigROIlist, channel, view);

  // step III.3: verify that the values were actually moved
  BOOST_CHECK(sigROIlistCopy.empty());

} // WireTestCustomConstructors()


//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is the name of the test; each step may have a
// number of checks and it will fail if any of them does.
//

BOOST_AUTO_TEST_CASE(WireDefaultConstructor) {
  WireTestDefaultConstructor();
}

BOOST_AUTO_TEST_CASE(WireCustomConstructors) {
  WireTestCustomConstructors();
}
