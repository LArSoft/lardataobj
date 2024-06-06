/**
 * @file    Hit_test.cc
 * @brief   Simple test on a recob::Hit object
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    20150113
 * @version 1.0
 *
 * This test simply creates recob::Hit objects and verifies that the values it
 * can access are the right ones.
 *
 * See http://www.boost.org/libs/test for the Boost test library home page.
 *
 * Timing:
 * version 1.0: ~1.5" (debug mode)
 */

// C/C++ standard library
#include <numeric> // std::accumulate

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
#define BOOST_TEST_MODULE (hit_test)
#include "boost/test/unit_test.hpp"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"  // raw::ChannelID_t, raw::TDCtick_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t
#include "lardataobj/RecoBase/Hit.h"

//------------------------------------------------------------------------------
//--- Test code
//

void CheckHit(recob::Hit const& hit,
              raw::ChannelID_t channel,
              raw::TDCtick_t start_tick,
              raw::TDCtick_t end_tick,
              float peak_time,
              float sigma_peak_time,
              float rms,
              float peak_amplitude,
              float sigma_peak_amplitude,
              float ROIsummedADC,
              float HitsummedADC,
              float hit_integral,
              float hit_sigma_integral,
              short int multiplicity,
              short int local_index,
              float goodness_of_fit,
              int dof,
              geo::View_t view,
              geo::SigType_t signal_type,
              geo::WireID wireID)
{

  // verify that the values are as expected
  // - channel ID
  BOOST_TEST(hit.Channel() == channel);

  // - view
  BOOST_TEST(hit.View() == view);

  // - signal type
  BOOST_TEST(hit.SignalType() == signal_type);

  // - start and end tick
  BOOST_TEST(hit.StartTick() == start_tick);
  BOOST_TEST(hit.EndTick() == end_tick);

  // - peak
  BOOST_TEST(hit.PeakTime() == peak_time);
  BOOST_TEST(hit.SigmaPeakTime() == sigma_peak_time);
  BOOST_TEST(hit.PeakAmplitude() == peak_amplitude);
  BOOST_TEST(hit.SigmaPeakAmplitude() == sigma_peak_amplitude);

  // the following comparisons are at 0.01%
  BOOST_CHECK_CLOSE(hit.PeakTimePlusRMS(), peak_time + rms, 0.01);
  BOOST_CHECK_CLOSE(hit.PeakTimeMinusRMS(), peak_time - rms, 0.01);

  for (float shift : {0.0, 0.5, 1.0, 1.5, 2.0, 2.2}) {

    const float time_up = peak_time + shift * rms;
    const float time_down = peak_time - shift * rms;

    if (time_up == 0.) {
      BOOST_CHECK_SMALL(hit.PeakTimePlusRMS(shift), 0.01F);
      BOOST_CHECK_SMALL(hit.PeakTimeMinusRMS(-shift), 0.01F);
    }
    else {
      BOOST_CHECK_CLOSE(hit.PeakTimePlusRMS(shift), time_up, 0.01F);
      BOOST_CHECK_CLOSE(hit.PeakTimeMinusRMS(-shift), time_up, 0.01F);
    }

    if (time_down == 0.) {
      BOOST_CHECK_SMALL(hit.PeakTimePlusRMS(shift), 0.01F);
      BOOST_CHECK_SMALL(hit.PeakTimeMinusRMS(-shift), 0.01F);
    }
    else {
      BOOST_CHECK_CLOSE(hit.PeakTimeMinusRMS(shift), time_down, 0.01F);
      BOOST_CHECK_CLOSE(hit.PeakTimePlusRMS(-shift), time_down, 0.01F);
    }

    if (rms > 0.) {
      if (shift == 0.) {
        BOOST_CHECK_SMALL(hit.TimeDistanceAsRMS(time_up), 0.01F);
        BOOST_CHECK_SMALL(hit.TimeDistanceAsRMS(time_down), 0.01F);
      }
      else {
        BOOST_CHECK_CLOSE(hit.TimeDistanceAsRMS(time_up), shift, 0.01F);
        BOOST_CHECK_CLOSE(hit.TimeDistanceAsRMS(time_down), -shift, 0.01F);
      }
    } // if rms is not 0

  } // for

  // - width
  BOOST_TEST(hit.RMS() == rms);

  // - charge
  BOOST_TEST(hit.ROISummedADC() == ROIsummedADC);
  BOOST_TEST(hit.Integral() == hit_integral);
  BOOST_TEST(hit.SigmaIntegral() == hit_sigma_integral);

  // - multiplicity
  BOOST_TEST(hit.Multiplicity() == multiplicity);
  BOOST_TEST(hit.LocalIndex() == local_index);
  BOOST_TEST(((hit.LocalIndex() < hit.Multiplicity()) || (hit.LocalIndex() == -1)));

  // - fit quality
  BOOST_TEST(hit.GoodnessOfFit() == goodness_of_fit);
  BOOST_TEST(hit.DegreesOfFreedom() == dof);

  // - wire ID
  BOOST_TEST(hit.WireID() == wireID);

} // CheckHit()

void HitTestDefaultConstructor()
{

  //
  // Part I: initialization of wire inputs
  //
  // these are the values expected for a default-constructed wire
  raw::ChannelID_t channel = raw::InvalidChannelID;
  raw::TDCtick_t start_tick = 0;
  raw::TDCtick_t end_tick = 0;
  float peak_time = 0.0;
  float sigma_peak_time = -1.0;
  float rms = 0.0;
  float peak_amplitude = 0.0;
  float sigma_peak_amplitude = -1.0;
  float ROIsummedADC = 0.0;
  float HitsummedADC = 0.0;
  float hit_integral = 0.0;
  float hit_sigma_integral = -1.0;
  short int multiplicity = 0;
  short int local_index = -1;
  float goodness_of_fit = 0.0;
  int dof = -1;
  geo::View_t view = geo::kUnknown;
  geo::SigType_t signal_type = geo::kMysteryType;
  geo::WireID wireID;

  //
  // Part II: default constructor
  //
  // step II.1: create a hit with the default constructor
  recob::Hit hit;

  // step II.2: verify that the values are as expected
  CheckHit(hit,
           channel,
           start_tick,
           end_tick,
           peak_time,
           sigma_peak_time,
           rms,
           peak_amplitude,
           sigma_peak_amplitude,
           ROIsummedADC,
           HitsummedADC,
           hit_integral,
           hit_sigma_integral,
           multiplicity,
           local_index,
           goodness_of_fit,
           dof,
           view,
           signal_type,
           wireID);

} // HitTestDefaultConstructor()

void HitTestCustomConstructors()
{

  //
  // Part I: initialization of wire inputs
  //
  // these are the values expected for a default-constructed wire
  raw::ChannelID_t channel = raw::InvalidChannelID;
  std::vector<float> signal({
    0.2,   1.2,   1.5,  2.0,  4.7,  12.8, 30.7, 85.2, 132.1, 154.7,
    147.4, 127.0, 86.7, 29.3, 11.8, 4.5,  2.2,  1.5,  1.0,   0.3,
  });
  raw::TDCtick_t start_tick = 730;
  raw::TDCtick_t end_tick = start_tick + signal.size();
  float peak_time = start_tick + signal.size() / 2 + 0.7;
  float sigma_peak_time = 1.3;
  float rms = signal.size() / 2.7;
  float peak_amplitude = signal[signal.size() / 2] - 0.3;
  float sigma_peak_amplitude = 2.3;
  float ROIsummedADC = std::accumulate(signal.begin(), signal.end(), 0.0);
  float HitsummedADC =
    std::accumulate(signal.begin(), signal.end(), 0.0); //NOT SURE WHAT IT SHOULD BE...
  float hit_integral = ROIsummedADC * 0.97;
  float hit_sigma_integral = peak_amplitude / 10.;
  short int multiplicity = 2;
  short int local_index = 1;
  float goodness_of_fit = 0.95;
  int dof = signal.size() - 3;
  geo::View_t view = geo::kU;
  geo::SigType_t signal_type = geo::kCollection;
  geo::WireID wireID(0, 1, 2, 546);

  //
  // Part II: complete constructor
  //
  // step II.1: create a hit with the signal-copying constructor
  recob::Hit hit1(channel,
                  start_tick,
                  end_tick,
                  peak_time,
                  sigma_peak_time,
                  rms,
                  peak_amplitude,
                  sigma_peak_amplitude,
                  ROIsummedADC,
                  HitsummedADC,
                  hit_integral,
                  hit_sigma_integral,
                  multiplicity,
                  local_index,
                  goodness_of_fit,
                  dof,
                  view,
                  signal_type,
                  wireID);

  // step II.2: verify that the values are as expected
  CheckHit(hit1,
           channel,
           start_tick,
           end_tick,
           peak_time,
           sigma_peak_time,
           rms,
           peak_amplitude,
           sigma_peak_amplitude,
           ROIsummedADC,
           HitsummedADC,
           hit_integral,
           hit_sigma_integral,
           multiplicity,
           local_index,
           goodness_of_fit,
           dof,
           view,
           signal_type,
           wireID);

} // HitTestCustomConstructors()

//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is the name of the test; each step may have a
// number of checks and it will fail if any of them does.
//

BOOST_AUTO_TEST_CASE(HitDefaultConstructor)
{
  HitTestDefaultConstructor();
}

BOOST_AUTO_TEST_CASE(HitCustomConstructors)
{
  HitTestCustomConstructors();
}
