/**
 * \file MCMiniPart.h
 *
 * \ingroup MCBase
 *
 * \brief Class def header for MCMiniPart data container
 *
 * @author Temigo
 */

/** \addtogroup MCBase

    @{*/
#ifndef MCMINIPART_H
#define MCMINIPART_H

// LArSoft
//#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCLimits.h" // kINVALID_X
//#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h" // simb::Origin_t

// ROOT libraries
#include <TVector3.h>

// STL
//#include <set>
#include <utility> // std::pair<>
#include <vector>

#include "TLorentzVector.h"

namespace sim
{

  class MCMiniPart {

  public:

    /// Default constructor
    MCMiniPart() {Reset();}

    /// Default destructor
    //virtual ~MCMiniPart() = default;

    void Reset() { ResetData(); }

    MCMiniPart(MCMiniPart const &) = default; // copy constructor
    MCMiniPart& operator=( const MCMiniPart &) = default;
    MCMiniPart(MCMiniPart&&) = default;
    MCMiniPart& operator= (MCMiniPart&&) = default;


    // Getters
    simb::Origin_t Origin () const { return _origin; }
    int            PdgCode() const { return _pdgcode; }
    unsigned int   TrackID() const { return _track_id; }
    const std::string&    Process() const { return _process; }
    unsigned int   Mother()  const { return _mother; }
    unsigned int   Ancestor() const { return _ancestor; }
    TLorentzVector StartVtx() const { return _start_vtx; }
    TLorentzVector StartMom() const { return _start_mom; }
    TLorentzVector EndVtx()   const { return _end_vtx; }
    TLorentzVector EndMom()   const { return _end_mom; }
    std::vector<std::pair<TLorentzVector,TLorentzVector> > DetPath() const { return _det_path; }
    std::vector<unsigned int> Daughters() const { return _daughters; }

    // Setters
    void Origin (simb::Origin_t o) { _origin = o; }
    void PdgCode(int id)           { _pdgcode = id; }
    void TrackID(unsigned int id)  { _track_id = id; }
    void Process(const std::string &name) { _process = name; }
    void Mother(unsigned int id)   { _mother = id; }
    void Ancestor(unsigned int id) { _ancestor = id; }
    void StartVtx(const TLorentzVector& vtx) { _start_vtx = vtx; }
    void StartMom(const TLorentzVector& mom) { _start_mom = mom; }
    void EndVtx(const TLorentzVector& vtx)   { _end_vtx = vtx; }
    void EndMom(const TLorentzVector& mom)   { _end_mom = mom; }
    void DetPath(const std::vector<std::pair<TLorentzVector,TLorentzVector> >& p) { _det_path = p; }
    void DetPath(std::vector<std::pair<TLorentzVector,TLorentzVector> >&& p) { _det_path = std::move(p); }
    void Daughters(const std::vector<unsigned int>& d) { _daughters = d; }
    void Daughters(std::vector<unsigned int>&& d) { _daughters = std::move(d); }

  protected:
    unsigned int   _track_id;
    std::string    _process;
    unsigned int   _mother;
    unsigned int   _ancestor;
    int            _pdgcode;
    TLorentzVector _start_vtx;
    TLorentzVector _start_mom;
    TLorentzVector _end_vtx;
    TLorentzVector _end_mom;
    std::vector<std::pair<TLorentzVector,TLorentzVector> > _det_path;
    std::vector<unsigned int> _daughters;
    ::simb::Origin_t _origin;

    void ResetData(){
      _track_id = _mother = _ancestor = kINVALID_UINT;
      _pdgcode  = kINVALID_INT;
      _process  = "";
      _origin   = ::simb::kUnknown;

      TLorentzVector invalid(kINVALID_DOUBLE,
           kINVALID_DOUBLE,
           kINVALID_DOUBLE,
           kINVALID_DOUBLE);
      _start_vtx = invalid;
      _start_mom = invalid;
      _end_vtx = invalid;
      _end_mom = invalid;
      _daughters.clear();
      _det_path.clear();
    }
  };
} // namespace sim

#endif
