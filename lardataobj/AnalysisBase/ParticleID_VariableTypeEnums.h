#ifndef PARTICLEID_VARIABLETYPEENUMS_H
#define PARTICLEID_VARIABLETYPEENUMS_H

namespace anab{

  enum kVariableType{
    kGOF,
    kLikelihood,
    kLogL,
    kScore,
    kPIDA,
    kdEdxtruncmean,
    kdQdxtruncmean,
    kTrackLength,
    kEdeposited,
    kEbyRange,
    kNotSet
  }; ///< Enum to define type of variable returned by a ParticleID algorithm

  enum kTrackDir{
    kForward,
    kBackward,
    kNoDirection
  }; ///< Enum to define track direction assumed by a ParticleID algorithm. Forward means the direction in which the track is reconstructed; backward means the opposite direction to the reconstruction.

}

#endif
