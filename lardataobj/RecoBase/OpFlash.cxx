////////////////////////////////////////////////////////////////////////
//
// \brief Definition of OpFlash reconstruction object
//
// \author bjpjones@mit.edu
//         cschiu@mit.edu
//
////////////////////////////////////////////////////////////////////////

#include "lardataobj/RecoBase/OpFlash.h"

#include <algorithm> // std::max()
#include <numeric>   // std::accumulate()
#include <ostream>
#include <utility> // std::move()

namespace recob {

  //----------------------------------------------------------------------
  OpFlash::OpFlash(double time,
                   double timewidth,
                   double abstime,
                   unsigned int frame,
                   std::vector<double> PEperOpDet,
                   bool InBeamFrame,
                   int onBeamTime,
                   double FastToTotal,
                   double xCenter,
                   double xWidth,
                   double yCenter,
                   double yWidth,
                   double zCenter,
                   double zWidth,
                   std::vector<double> WireCenters,
                   std::vector<double> WireWidths)
    : fTime{time}
    , fTimeWidth{timewidth}
    , fAbsTime{abstime}
    , fFrame{frame}
    , fPEperOpDet{std::move(PEperOpDet)}
    , fWireCenters{std::move(WireCenters)}
    , fWireWidths{std::move(WireWidths)}
    , fXCenter{xCenter}
    , fXWidth{xWidth}
    , fYCenter{yCenter}
    , fYWidth{yWidth}
    , fZCenter{zCenter}
    , fZWidth{zWidth}
    , fFastToTotal{FastToTotal}
    , fInBeamFrame{InBeamFrame}
    , fOnBeamTime{onBeamTime}
  {}

  //----------------------------------------------------------------------
  OpFlash::OpFlash(double time,
                   double timewidth,
                   double abstime,
                   unsigned int frame,
                   std::vector<double> PEperOpDet,
                   bool InBeamFrame,
                   int onBeamTime,
                   double FastToTotal,
                   double yCenter,
                   double yWidth,
                   double zCenter,
                   double zWidth,
                   std::vector<double> WireCenters,
                   std::vector<double> WireWidths)
    : OpFlash{time,
              timewidth,
              abstime,
              frame,
              std::move(PEperOpDet),
              InBeamFrame,
              onBeamTime,
              FastToTotal,
              NoCenter,
              NoCenter,
              yCenter,
              yWidth,
              zCenter,
              zWidth,
              std::move(WireCenters),
              std::move(WireWidths)}
  {}

  //----------------------------------------------------------------------
  bool operator<(const OpFlash& a, const OpFlash& b)
  {
    return a.Time() < b.Time();
  }

  //----------------------------------------------------------------------
  double OpFlash::TotalPE() const
  {
    return std::accumulate(fPEperOpDet.begin(), fPEperOpDet.end(), 0.0);
  }

  //----------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& out, OpFlash const& flash)
  {

    out << "Optical flash at " << flash.Time() << " us (absolute: " << flash.AbsTime() << " us) "
        << flash.TimeWidth() << " us long, centered at (";
    if (flash.hasXCenter())
      out << flash.XCenter();
    else
      out << "?";
    out << ", " << flash.YCenter() << ", " << flash.ZCenter() << ") cm of extension +/- (";
    if (flash.hasXCenter())
      out << flash.XWidth();
    else
      out << "?";
    out << ", " << flash.YWidth() << ", " << flash.ZWidth() << ") cm ";
    if (!flash.hasXCenter()) out << " (x: n/a)";

    out << "\n  " << (flash.InBeamFrame() ? "" : "not ")
        << "in beam frame (on beam time: " << flash.OnBeamTime() << ")";

    out << "\n  total p.e.: " << flash.TotalPE() << " (fast/total fraction: " << flash.FastToTotal()
        << ")";
    std::size_t const nPEs = flash.PEs().size();
    if (nPEs > 0) {
      out << " across " << nPEs << " channels; non-empty:";
      constexpr unsigned int pageSize = 8;
      unsigned int line = 0;
      for (std::size_t i = 0; i < nPEs; ++i) {
        double const pe = flash.PE(i);
        if (pe == 0.0) continue;
        if (line-- == 0) {
          line = pageSize;
          out << "\n  ";
        }
        out << "  [" << i << "] " << pe;
      }
    }
    else {
      out << ", no channels recorded;";
    }

    std::size_t const nWireCenters = flash.WireCenters().size();
    std::size_t const nWireWidths = flash.WireWidths().size();
    std::size_t const nWires = std::max(nWireCenters, nWireWidths);
    out << "\n  ";
    if (nWires) {
      out << nWires << " wires recorded (center +/- width):";
      constexpr unsigned int pageSize = 5;
      unsigned int line = 0;
      for (std::size_t i = 0; i < nWires; ++i) {
        if (line-- == 0) {
          line = pageSize;
          out << "\n  ";
        }
        out << "  [" << i << "] (";
        if (i < nWireCenters)
          out << flash.WireCenters()[i];
        else
          out << "n/a";
        out << " +/- ";
        if (i < nWireWidths)
          out << flash.WireWidths()[i];
        else
          out << "n/a";
        out << ")";
      }
    }
    else
      out << "no wire information recorded";

    return out;
  } // std::ostream& operator<<(std::ostream&, OpFlash const&)

  //----------------------------------------------------------------------

}
