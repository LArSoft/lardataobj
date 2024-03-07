////////////////////////////////////////////////////////////////////////
//
// A timestamp data product to be associated with raw::RawDigit
// Tom Junk, June 26, 2018
//
////////////////////////////////////////////////////////////////////////

#ifndef RDTIMESTAMPDATAPRODUCT_H
#define RDTIMESTAMPDATAPRODUCT_H

// uint64_t is typedef unsigned long long on macOS but unsigned long on Linux
// We use the root type ULong64_t from RtypesCore.h instead
#include "RtypesCore.h"
#include <stdint.h>

namespace raw {

  class RDTimeStamp {

  public:
    RDTimeStamp(ULong64_t tstamp = 0, uint16_t flags = 0)
    {
      fTimeStamp = tstamp;
      fFlags = flags;
    } // Constructor

    void SetTimeStamp(ULong64_t tstamp);
    void SetTimeStamp_Low(uint32_t tslow);
    void SetTimeStamp_High(uint32_t tshigh);
    ULong64_t GetTimeStamp() const;
    uint32_t GetTimeStamp_Low() const;
    uint32_t GetTimeStamp_High() const;
    void SetFlags(uint16_t flags);
    uint16_t GetFlags() const;

  private:
    ULong64_t fTimeStamp;
    uint16_t fFlags;
  };

  inline void RDTimeStamp::SetTimeStamp(ULong64_t tstamp)
  {
    fTimeStamp = tstamp;
  }
  inline void RDTimeStamp::SetTimeStamp_Low(uint32_t tslow)
  {
    fTimeStamp &= 0xFFFFFFFF00000000;
    fTimeStamp |= tslow;
  }
  inline void RDTimeStamp::SetTimeStamp_High(uint32_t tshigh)
  {
    fTimeStamp &= 0xFFFFFFFF;
    fTimeStamp |= ((ULong64_t)tshigh << 32);
  }
  inline ULong64_t RDTimeStamp::GetTimeStamp() const
  {
    return fTimeStamp;
  }
  inline uint32_t RDTimeStamp::GetTimeStamp_Low() const
  {
    return fTimeStamp & 0xFFFFFFFF;
  }
  inline uint32_t RDTimeStamp::GetTimeStamp_High() const
  {
    return fTimeStamp >> 32;
  }
  inline void RDTimeStamp::SetFlags(uint16_t flags)
  {
    fFlags = flags;
  }
  inline uint16_t RDTimeStamp::GetFlags() const
  {
    return fFlags;
  }

}

#endif // RDTIMESTAMPDATAPRODUCT_H
