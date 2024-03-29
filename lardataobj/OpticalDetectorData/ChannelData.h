// OpticalDetectorData/ChannelData.h
// William Seligman <seligman@nevis.columbia.edu>

// The ADC counts associated with a particular channel.

#ifndef OpticalDetectorData_ChannelData_h
#define OpticalDetectorData_ChannelData_h

// LArSoft includes
#include "lardataobj/OpticalDetectorData/OpticalTypes.h"

// C++ includes
#include <functional> // so we can redefine less<> below
#include <limits>
#include <vector>

namespace optdata {

  class ChannelData : public std::vector<ADC_Count_t> {
  public:
    // Simple constructors/destructors.
    // Just in case the user forgets to supply the default channel, use
    // a garbage value to indicate that there's a problem.
    // To save on memory reallocations, offer an option to specify the
    // the initial memory allocation of the channel vector.
    ChannelData(Channel_t chan = std::numeric_limits<Channel_t>::max(), size_type len = 0)
      : fm_optDetChannel(chan)
    {
      this->reserve(len);
    };

    ~ChannelData(){};

    // No "setter" for the channel number; you have to assign it when
    // you create a ChannelData object.
    Channel_t ChannelNumber() const { return fm_optDetChannel; }

  private:
    unsigned int fm_optDetChannel;
  };

  // In case we want to sort a collection of ChannelDatas (e.g.,
  // std::set<ChannelData>), here's the definition of the less-than
  // operator.
  bool operator<(const ChannelData& lhs, const ChannelData& rhs)
  {
    // Sort by channel.
    if (lhs.ChannelNumber() < rhs.ChannelNumber()) return true;
    return false;
  }

} // namespace optdata

// For no extra charge, include how to sort ChannelData*, just in
// case we want (for example) a std::set<ChannelData*>.
namespace std {
  template <>
  class less<optdata::ChannelData*> {
  public:
    bool operator()(const optdata::ChannelData* lhs, const optdata::ChannelData* rhs)
    {
      return (*lhs) < (*rhs);
    }
  };
} // std

#endif // OpticalDetectorData_ChannelData_h
