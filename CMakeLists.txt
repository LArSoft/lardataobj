# ======================================================================
#  lardataobj main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/lardataobj/ups/setup_for_development <-d|-p>
#  buildtool [-I /install/path] -bt[i] -jN
# ======================================================================

cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(lardataobj LANGUAGES CXX)

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

include(CetCMakeEnv)
set(${PROJECT}_CMAKE_PROJECT_VERSION_STRING 09.05.00-rc0)
cet_cmake_env()

cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

find_package(cetlib_except REQUIRED PUBLIC)
find_package(cetlib REQUIRED PUBLIC)
find_package(messagefacility REQUIRED PUBLIC)
find_package(canvas REQUIRED PUBLIC)
find_package(canvas_root_io REQUIRED)
find_package(nusimdata REQUIRED PUBLIC)
find_package(Boost REQUIRED)
find_package(ROOT COMPONENTS Core Physics GenVector Matrix REQUIRED PUBLIC)
# these are minimum required versions, not the actual product versions
find_ups_product( larcoreobj )
find_ups_product( larcorealg )

# macros for artdaq_dictionary and simple_plugin
include(ArtDictionary)
include(CetMake)
include(BasicPlugin)

# source
add_subdirectory(lardataobj)

# tests
add_subdirectory(test)

# ups - table and config files
add_subdirectory(ups)

# packaging utility
include(UseCPack)
