////////////////////////////////////////////////////////////////////////
//
// OpDetPulse class
//
// Ben Jones
//  bjpjones@mit.edu
//
// Christie Chiu
//  cschiu@mit.edu
//
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RawData/OpDetPulse.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace raw {

  //----------------------------------------------------------------------
  OpDetPulse::OpDetPulse()
  {
    /*    mf::LogWarning("OpDetPulse") <<"Warning : OpDetPulse default constructor called without "
				 <<"specification of OpChannel.  This is not recommended";
    */
  }

  //----------------------------------------------------------------------
  OpDetPulse::OpDetPulse(unsigned short OpChannel,
                         std::vector<short> Waveform,
                         unsigned int PMTFrame,
                         unsigned int FirstSample)
  {
    fWaveform = Waveform;
    fOpChannel = OpChannel;
    fPMTFrame = PMTFrame;
    fFirstSample = FirstSample;
  }

  //----------------------------------------------------------------------
  OpDetPulse::OpDetPulse(unsigned short OpChannel)
  {
    fOpChannel = OpChannel;
  }

  //----------------------------------------------------------------------
  int OpDetPulse::Integral(unsigned short startbin, unsigned short endbin)
  {
    // endbin set to zero or not supplied means go to end of vector
    if (endbin == 0) endbin = Samples();

    // Some error trapping for stupid integration limits
    if (endbin < startbin)
      throw cet::exception("OpDetPulse")
        << " Trying to integrate between illegal limits : " << startbin << " " << endbin << "\n";
    if (endbin > Samples())
      throw cet::exception("OpDetPulse")
        << " Trying to integrate between illegal limits : " << startbin << " " << endbin << "\n";

    unsigned int TotalIntegral = 0;
    // Loop over vector to return the integral
    for (size_t i = startbin; i != endbin; ++i) {
      TotalIntegral += fWaveform.at(i);
    }

    return TotalIntegral;
  }

}
