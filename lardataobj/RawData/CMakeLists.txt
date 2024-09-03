cet_make_library(SOURCE
  AuxDetDigit.cxx
  BeamInfo.cxx
  DAQHeader.cxx
  ExternalTrigger.cxx
  OpDetPulse.cxx
  raw.cxx
  RawDigit.cxx
  RDTimeStamp.h
  TriggerData.cxx
  LIBRARIES
  PUBLIC
  Boost::boost
  ROOT::Core
  larcoreobj::headers
  PRIVATE
  messagefacility::MF_MessageLogger
  cetlib_except::cetlib_except
)

build_dictionary(DICTIONARY_LIBRARIES
  lardataobj::RawData
  canvas::canvas
)

install_headers()
install_source()
