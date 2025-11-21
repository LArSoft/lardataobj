cet_make_library(LIBRARY_NAME TrackingTypes INTERFACE
  SOURCE TrackingTypes.h
  LIBRARIES INTERFACE
  larcoreobj::geo_vectors
  ROOT::GenVector
  ROOT::Matrix
  ROOT::Physics
)

cet_make_library(LIBRARY_NAME TrackFitHitInfo INTERFACE
  SOURCE TrackFitHitInfo.h
  LIBRARIES INTERFACE
  lardataobj::TrackingTypes
  larcoreobj::SimpleTypesAndConstants
)

cet_make_library(SOURCE
  Cluster.cxx
  Edge.cxx
  EndPoint2D.cxx
  Event.cxx
  Hit.cxx
  OpFlash.cxx
  OpHit.cxx
  PCAxis.cxx
  PFParticle.cxx
  PFParticleMetadata.cxx
  Seed.cxx
  Shower.cxx
  Slice.cxx
  SpacePoint.cxx
  Track.cxx
  TrackingPlane.cxx
  TrackTrajectory.cxx
  Trajectory.cxx
  TrajectoryPointFlags.cxx
  Vertex.cxx
  Wire.cxx
  OpWaveform.cxx
  LIBRARIES
  PUBLIC
  lardataobj::TrackingTypes
  larcoreobj::SimpleTypesAndConstants
  larcoreobj::geo_vectors
  ROOT::Core
  ROOT::Physics
  PRIVATE
  messagefacility::MF_MessageLogger
)

build_dictionary(DICTIONARY_LIBRARIES
  lardataobj::RecoBase
  lardataobj::TrackFitHitInfo
  lardataobj::TrackingTypes
  canvas::canvas
)

add_subdirectory(AssnsDicts)
add_subdirectory(TrackingDicts)

install_headers()
install_source(EXTRAS RecoBase.dox)
