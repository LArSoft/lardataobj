/// OpticalRawDigit.h
/// Kazuhiro Terao <kazuhiro@nevis.columbia.edu>
///
/// Equivalent of raw::RawDigit (waveform for TPC) for Optical system
/// except the code is a complete copy from FIFOChannel written by William Seligman.
/// William Seligman <seligman@nevis.columbia.edu>
///

#ifndef optdata_OpticalRawDigit_h
#define optdata_OpticalRawDigit_h

// LArSoft includes
#include "lardataobj/OpticalDetectorData/FIFOChannel.h"

// C++ includes
#include <functional> // so we can redefine less<> below
#include <limits>
#include <vector>

namespace optdata {

  class OpticalRawDigit : public FIFOChannel {
  public:
    // Simple constructor/destructor.
    OpticalRawDigit(Optical_Category_t category = kUndefined,
                    TimeSlice_t time = 0,
                    Frame_t frame = 0,
                    Channel_t channel = std::numeric_limits<Channel_t>::max(),
                    size_type len = 0)
      : FIFOChannel(category, time, frame, channel, len){};

    ~OpticalRawDigit(){};
  };

  // In case we want to sort a collection of OpticalRawDigits (e.g.,
  // std::set<OpticalRawDigit>), here's the definition of the less-than
  // operator.
  bool operator<(const OpticalRawDigit& lhs, const OpticalRawDigit& rhs)
  {
    // Sort by channel, frame number, and time associated with the first bin.
    if (lhs.ChannelNumber() < rhs.ChannelNumber() && lhs.Frame() < rhs.Frame() &&
        lhs.TimeSlice() < rhs.TimeSlice())
      return true;
    return false;
  }

} // namespace optdata

// For no extra charge, include how to sort OpticalRawDigit*, just in
// case we want (for example) a std::set<OpticalRawDigit*>.
namespace std {
  template <>
  class less<optdata::OpticalRawDigit*> {
  public:
    bool operator()(const optdata::OpticalRawDigit* lhs, const optdata::OpticalRawDigit* rhs)
    {
      return (*lhs) < (*rhs);
    }
  };
} // std

#endif // optdata_OpticalRawDigit_h
