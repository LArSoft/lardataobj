/**
 * @file    SimChannel_test.cc
 * @brief   Limited test on a `sim::SimChannel` object
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    December 25, 2025
 */

// LArSoft libraries
#include "lardataobj/Simulation/SimChannel.h"

// C++ standard library
#include <array>
#include <initializer_list>
#include <random>

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
#define BOOST_TEST_MODULE (simchannel_test)
#include "boost/test/unit_test.hpp"

//------------------------------------------------------------------------------
//--- Test code
//
/**
 * Returns a `sim::SimChannel` with data in specified TDC intervals.
 *
 * On each interval, `{ startTick, endTick, nIDE }`, `nIDE` IDE are added
 * randomly to TDC between `startTick` (included) and `endTick` (excluded).
 * All IDE are added identical (final track 1000 at origin, 1 MeV and 500
 * electrons) but all different original tracks.
 */
sim::SimChannel makeSimChannelWithDataBetween(
  std::initializer_list<std::array<unsigned int, 3>> intervals)
{
  sim::SimChannel sc{12345};

  std::default_random_engine rndEng{123456789}; // make this reproducible

  auto addIDEat = [&sc](int trackID, unsigned int tick) {
    static std::array const pos = {0., 0., 0.};
    sc.AddIonizationElectrons(1000, tick, 500., pos.data(), 1.0, trackID);
  };

  for (auto const [startTick, endTick, nIDE] : intervals) {
    if (startTick >= endTick) continue;
    std::uniform_int_distribution<> extractTick(startTick, endTick - 1);
    for (unsigned int ID = 0; ID < nIDE; ++ID)
      addIDEat(ID, extractTick(rndEng));
  } // for intervals

  return sc;
} // makeSimChannelWithDataBetween()

// -----------------------------------------------------------------------------
template <typename Iter>
void CheckIDE(Iter itIDE, Iter begin, Iter end, unsigned int tick)
{

  Iter it = begin;
  for (; it != itIDE; ++it)
    BOOST_TEST(it->first < tick);
  for (; it != end; ++it)
    BOOST_TEST(it->first >= tick);
}

template <typename Iter>
void CheckIDEspan(Iter beginSpan,
                  Iter endSpan,
                  Iter begin,
                  Iter end,
                  unsigned int beginTick,
                  unsigned int endTick)
{

  Iter it = begin;
  for (; it != beginSpan; ++it)
    BOOST_TEST(it->first < beginTick);
  for (; it != endSpan; ++it) {
    BOOST_TEST(it->first >= beginTick);
    BOOST_TEST(it->first < endTick);
  }
  for (; it != end; ++it)
    BOOST_TEST(it->first >= endTick);
}

// -----------------------------------------------------------------------------
void SimChannel_findIDEatTime_test()
{

  /*
   * The promise:
   *
   * double chargeBetween(sim::SimChannel const& sc, int start, int ticks) {
   *   double charge = 0.0;
   *   for (sim::TDCIDE const& TDCIDE: sc.IDEsBetween(start, start + ticks)) {
   *     for (sim::IDE const& ide: TDCIDE.second)
   *       charge += ide.numElectrons;
   *   }
   *   return charge;
   * }
   *
   */

  sim::SimChannel const sc = makeSimChannelWithDataBetween({{1500, 1510, 15}, {1550, 1560, 10}});

  sim::TDCIDE const* firstIDE = &(sc.TDCIDEMap().front());
  sim::TDCIDE const* endIDE = firstIDE + sc.TDCIDEMap().size();

  BOOST_TEST(sc.findIDEatTime(1400) == firstIDE);
  BOOST_TEST(sc.findIDEatTime(1500) == firstIDE);
  sim::TDCIDE const* IDE1550 = sc.findIDEatTime(1550);
  BOOST_CHECK(IDE1550);
  if (IDE1550) {
    BOOST_TEST_INFO("SimChannel_findIDEatTime_test: findIDEatTime(1550)");
    CheckIDE(IDE1550, firstIDE, endIDE, 1550);
  }
  BOOST_TEST(sc.findIDEatTime(1510) == IDE1550);
  BOOST_TEST(sc.findIDEatTime(1520) == IDE1550);
  BOOST_TEST(sc.findIDEatTime(1530) == IDE1550);
  BOOST_CHECK(sc.findIDEatTime(std::prev(endIDE)->first));
  BOOST_CHECK(!sc.findIDEatTime(1560));
  BOOST_CHECK(!sc.findIDEatTime(1600));

} // SimChannel_findIDEatTime_test()

// -----------------------------------------------------------------------------
void SimChannel_IDEsBetween_test()
{

  sim::SimChannel const sc = makeSimChannelWithDataBetween({{1500, 1510, 15}, {1550, 1560, 10}});

  auto const checkIDEsBetween = [&sc](unsigned int start, unsigned int stop, unsigned int n) {
    BOOST_TEST_INFO("SimChannel_findIDEatTime_test: checkIDEs(" << start << ", " << stop << ", "
                                                                << n << ")");
    // test range
    auto const IDEspan = sc.IDEsBetween(start, stop);
    auto const& IDEmap = sc.TDCIDEMap();
    using std::begin, std::end;
    CheckIDEspan(begin(IDEspan), end(IDEspan), IDEmap.begin(), IDEmap.end(), start, stop);
    // test count
    std::size_t count = 0;
    for (sim::TDCIDE const& IDEs : IDEspan)
      count += IDEs.second.size();
    BOOST_TEST(count == n);
  };

  checkIDEsBetween(1400, 1450, 0);
  checkIDEsBetween(1400, 1500, 0);
  checkIDEsBetween(1500, 1510, 15);
  checkIDEsBetween(1400, 1510, 15);
  checkIDEsBetween(1400, 1520, 15);
  checkIDEsBetween(1400, 1550, 15);
  checkIDEsBetween(1510, 1540, 0);
  checkIDEsBetween(1510, 1550, 0);
  checkIDEsBetween(1520, 1540, 0);
  checkIDEsBetween(1520, 1550, 0);
  checkIDEsBetween(1550, 1560, 10);
  checkIDEsBetween(1510, 1560, 10);
  checkIDEsBetween(1520, 1560, 10);
  checkIDEsBetween(1520, 1580, 10);
  checkIDEsBetween(1550, 1580, 10);
  checkIDEsBetween(1560, 1580, 0);
  checkIDEsBetween(1570, 1580, 0);
  checkIDEsBetween(1500, 1560, 25);
  checkIDEsBetween(1500, 1680, 25);
  checkIDEsBetween(1400, 1560, 25);
  checkIDEsBetween(1400, 1580, 25);

} // SimChannel_IDEsBetween_test()

// -----------------------------------------------------------------------------
sim::SimChannel makeSimChannel()
{

  sim::SimChannel sc{1234};
  std::array<double, 3U> pos;

  // AddIonizationElectrons(track ID, TDC, nElectrons, location, energy, [origTrackID]);

  // * TDC 1500: two tracks
  pos = {1.0, -1.0, 2.0};
  sc.AddIonizationElectrons(10, 1500, 500.0, pos.data(), 0.02);
  sc.AddIonizationElectrons(10, 1500, 1000.0, pos.data(), 0.04);
  pos = {1.0, 5.0, 20.0};
  sc.AddIonizationElectrons(20, 1500, 1000.0, pos.data(), 0.04);

  // * TDC 1502: one track
  pos = {1.3, -1.0, 2.0};
  sc.AddIonizationElectrons(10, 1502, 2000.0, pos.data(), 0.08);

  return sc;
} // makeSimChannel()

double chargeBetween(sim::SimChannel const& sc, int start, int ticks)
{
  double charge = 0.0;
  for (sim::TDCIDE const& TDCIDE : sc.IDEsBetween(start, start + ticks)) {
    for (sim::IDE const& ide : TDCIDE.second)
      charge += ide.numElectrons;
  }
  return charge;
}

void SimChannel_IDEsBetween_documentation_test()
{

  /*
   * The promise:
   *
   * double chargeBetween(sim::SimChannel const& sc, int start, int ticks) {
   *   double charge = 0.0;
   *   for (sim::TDCIDE const& TDCIDE: sc.IDEsBetween(start, start + ticks)) {
   *     for (sim::IDE const& ide: TDCIDE.second)
   *       charge += ide.numElectrons;
   *   }
   *   return charge;
   * }
   *
   */

  sim::SimChannel const sc = makeSimChannel();

  BOOST_REQUIRE(sc.TDCIDEMap().front().first == 1500);
  BOOST_REQUIRE(sc.TDCIDEMap().back().first < 1600);
  BOOST_TEST(chargeBetween(sc, 1400, 1600) == 4500.0, 0.01 % boost::test_tools::tolerance());
  BOOST_TEST(chargeBetween(sc, 1500, 1600) == 4500.0, 0.01 % boost::test_tools::tolerance());
  BOOST_TEST(chargeBetween(sc, 1501, 1600) == 2000.0, 0.01 % boost::test_tools::tolerance());
  BOOST_TEST(chargeBetween(sc, 1502, 1600) == 2000.0, 0.01 % boost::test_tools::tolerance());
  BOOST_TEST(chargeBetween(sc, 1503, 1600) == 0.0);

} // SimChannel_IDEsBetween_documentation_test()

//------------------------------------------------------------------------------
//--- registration of tests
//
// Boost needs now to know which tests we want to run.
// Tests are "automatically" registered, hence the BOOST_AUTO_TEST_CASE()
// macro name. The argument is the name of the test; each step may have a
// number of checks and it will fail if any of them does.
//

BOOST_AUTO_TEST_CASE(SimChannel_lookup_testcase)
{
  SimChannel_findIDEatTime_test();
  SimChannel_IDEsBetween_test();
}

BOOST_AUTO_TEST_CASE(SimChannel_documentation_testcase)
{
  SimChannel_IDEsBetween_documentation_test();
}
