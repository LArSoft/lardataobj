////////////////////////////////////////////////////////////////////////////
//
// \brief 3D seed object for kalman tracking package
//                and bezier tracking algorithm
//
// \author Ben Jones, MIT
//         bjpjones@mit.edu
//
////////////////////////////////////////////////////////////////////////////
#include "lardataobj/RecoBase/Seed.h"
#include "TVector3.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include <iostream>
#include <stdlib.h>

namespace recob {

  //----------------------------------------------------------------------------
  Seed::Seed()
  {
    fIsValid = false;
  }

  //----------------------------------------------------------------------------
  Seed::Seed(double* Pt, double* Dir)
  {
    for (int i = 0; i != 3; i++) {
      fSeedPoint[i] = Pt[i];
      fSeedDirection[i] = Dir[i];
      fSeedPointError[i] = 0;
      fSeedDirectionError[i] = 0;
    }
    fIsValid = true;
  }

  //----------------------------------------------------------------------------
  Seed::Seed(double* Pt, double* Dir, double* PtErr, double* DirErr)
  {
    for (int i = 0; i != 3; i++) {
      fSeedPoint[i] = Pt[i];
      fSeedDirection[i] = Dir[i];
      fSeedPointError[i] = PtErr[i];
      fSeedDirectionError[i] = DirErr[i];
    }
    fIsValid = true;
  }

  //----------------------------------------------------------------------------
  void Seed::Print() const
  {
    std::cout << "Printing seed contents : " << fSeedPoint[0] << " " << fSeedPoint[1] << " "
              << fSeedPoint[2] << ", " << fSeedDirection[0] << " " << fSeedDirection[1] << " "
              << fSeedDirection[2] << std::endl;
  }

  //----------------------------------------------------------------------------
  bool Seed::IsValid() const
  {
    return fIsValid;
  }

  //----------------------------------------------------------------------------
  Seed Seed::Reverse()
  {
    double NewSeedDir[3];
    for (size_t n = 0; n != 3; ++n) {
      NewSeedDir[n] = -fSeedDirection[n];
    }
    return Seed(fSeedPoint, NewSeedDir, fSeedPointError, fSeedDirectionError);
  }

  //----------------------------------------------------------------------------
  void Seed::SetValidity(bool Validity)
  {
    fIsValid = Validity;
  }

  //----------------------------------------------------------------------------
  void Seed::GetDirection(double* rDir, double* rDirErr) const
  {
    for (int i = 0; i != 3; i++) {
      rDir[i] = fSeedDirection[i];
      if (rDirErr) rDirErr[i] = fSeedDirectionError[i];
    }
  }

  //----------------------------------------------------------------------------
  void Seed::GetPoint(double* rPt, double* rPtErr) const
  {
    for (int i = 0; i != 3; i++) {
      rPt[i] = fSeedPoint[i];
      if (rPtErr) rPtErr[i] = fSeedPointError[i];
    }
  }

  //----------------------------------------------------------------------------
  void Seed::SetDirection(double* Dir)
  {
    double Empty[3] = {0, 0, 0};
    SetDirection(Dir, Empty);
  }

  //----------------------------------------------------------------------------
  void Seed::SetPoint(double* Pt)
  {
    double Empty[3] = {0, 0, 0};
    SetPoint(Pt, Empty);
  }

  //----------------------------------------------------------------------------
  void Seed::SetDirection(double* Dir, double* Err)
  {
    for (int i = 0; i != 3; i++) {
      fSeedDirection[i] = Dir[i];
      fSeedDirectionError[i] = Err[i];
    }
    fIsValid = true;
  }

  //----------------------------------------------------------------------------
  void Seed::SetPoint(double* Pt, double* Err)
  {
    for (int i = 0; i != 3; i++) {
      fSeedPoint[i] = Pt[i];
      fSeedPointError[i] = Err[i];
    }
    fIsValid = true;
  }

  //----------------------------------------------------------------------------
  double Seed::GetLength() const
  {
    return pow(pow(fSeedDirection[0], 2) + pow(fSeedDirection[1], 2) + pow(fSeedDirection[2], 2),
               0.5);
  }

  //----------------------------------------------------------------------------
  double Seed::GetProjAngleDiscrepancy(Seed const& AnotherSeed) const
  {
    double OtherPt[3];
    double OtherPtErr[3];

    AnotherSeed.GetPoint(OtherPt, OtherPtErr);

    TVector3 OtherPtV(OtherPt[0], OtherPt[1], OtherPt[2]);
    TVector3 ThisDirV(fSeedDirection[0], fSeedDirection[1], fSeedDirection[2]);
    TVector3 ThisPtV(fSeedPoint[0], fSeedPoint[1], fSeedPoint[2]);

    return (OtherPtV - ThisPtV).Angle(ThisDirV.Unit());
  }

  //----------------------------------------------------------------------------
  void Seed::GetVectorBetween(Seed const& AnotherSeed, double* xyz) const
  {
    double xyzother[3], err[3];
    AnotherSeed.GetPoint(xyzother, err);

    xyz[0] = xyzother[0] - fSeedPoint[0];
    xyz[1] = xyzother[1] - fSeedPoint[1];
    xyz[2] = xyzother[2] - fSeedPoint[2];
  }

  //----------------------------------------------------------------------------
  double Seed::GetAngle(Seed const& AnotherSeed) const
  {
    double OtherDir[3];
    double OtherDirErr[3];
    AnotherSeed.GetDirection(OtherDir, OtherDirErr);

    double OtherMag = AnotherSeed.GetLength();
    double ThisMag = GetLength();

    /*
    std::cout<<"Seed angle calc: " <<
      OtherMag<< " " << ThisMag << " " <<
      ( (OtherDir[0]*fSeedDirection[0] +
   OtherDir[1]*fSeedDirection[1] +
   OtherDir[2]*fSeedDirection[2] ) / (ThisMag*OtherMag)) <<
      std::endl;
    */

    // Need a tiny offset, as acos(1.0) gives unpredictable results due to
    // floating point precision
    double eta = 0.00001;

    return std::acos((OtherDir[0] * fSeedDirection[0] + OtherDir[1] * fSeedDirection[1] +
                      OtherDir[2] * fSeedDirection[2] - eta) /
                     (ThisMag * OtherMag));
  }

  //----------------------------------------------------------------------------
  double Seed::GetProjDiscrepancy(Seed const& AnotherSeed) const
  {
    double OtherPt[3];
    double OtherPtErr[3];

    AnotherSeed.GetPoint(OtherPt, OtherPtErr);

    TVector3 OtherPtV(OtherPt[0], OtherPt[1], OtherPt[2]);
    TVector3 ThisDirV(fSeedDirection[0], fSeedDirection[1], fSeedDirection[2]);
    TVector3 ThisPtV(fSeedPoint[0], fSeedPoint[1], fSeedPoint[2]);

    return ((OtherPtV - ThisPtV) - ThisDirV.Unit() * (ThisDirV.Unit().Dot(OtherPtV - ThisPtV)))
      .Mag();
  }

  //----------------------------------------------------------------------------
  double Seed::GetDistance(Seed const& AnotherSeed) const
  {
    double OtherPt[3];
    double OtherPtErr[3];

    AnotherSeed.GetPoint(OtherPt, OtherPtErr);

    return pow(pow(OtherPt[0] - fSeedPoint[0], 2) + pow(OtherPt[1] - fSeedPoint[1], 2) +
                 pow(OtherPt[2] - fSeedPoint[2], 2),
               0.5);
  }

  //----------------------------------------------------------------------------
  double Seed::GetDistanceFrom(recob::SpacePoint const& SomePoint) const
  {
    double SPxyz[3];
    SPxyz[0] = SomePoint.XYZ()[0];
    SPxyz[1] = SomePoint.XYZ()[1];
    SPxyz[2] = SomePoint.XYZ()[2];

    //  std::cout<<"Seed Dir Vec " << fSeedDirection[0]<<" " <<fSeedDirection[1]<< " " << fSeedDirection[2]<<std::endl;

    double ThisSeedLength =
      pow(pow(fSeedDirection[0], 2) + pow(fSeedDirection[1], 2) + pow(fSeedDirection[2], 2), 0.5);

    double SPProjOnSeed = (fSeedDirection[0] * (SPxyz[0] - fSeedPoint[0]) +
                           fSeedDirection[1] * (SPxyz[1] - fSeedPoint[1]) +
                           fSeedDirection[2] * (SPxyz[2] - fSeedPoint[2])) /
                          ThisSeedLength;

    // std::cout<<"proj : " <<SPProjOnSeed<<std::endl;
    // std::cout<<"Seed len :" << ThisSeedLength<<std::endl;

    if (SPProjOnSeed > (ThisSeedLength)) {
      //  std::cout<<"Seed over end"<<std::endl;
      return pow(pow(fSeedPoint[0] + fSeedDirection[0] - SPxyz[0], 2) +
                   pow(fSeedPoint[1] + fSeedDirection[1] - SPxyz[1], 2) +
                   pow(fSeedPoint[2] + fSeedDirection[2] - SPxyz[2], 2),
                 0.5);
    }
    else if (SPProjOnSeed < (0 - ThisSeedLength)) {
      //  std::cout<<"Seed under end"<<std::endl;
      return pow(pow(fSeedPoint[0] - fSeedDirection[0] - SPxyz[0], 2) +
                   pow(fSeedPoint[1] - fSeedDirection[1] - SPxyz[1], 2) +
                   pow(fSeedPoint[2] - fSeedDirection[2] - SPxyz[2], 2),
                 0.5);
    }
    else {
      //  std::cout<<"Seed valid region"<<std::endl;
      double crossprod[3];
      CrossProd(fSeedPoint[0] + fSeedDirection[0] - SPxyz[0],
                fSeedPoint[1] + fSeedDirection[1] - SPxyz[1],
                fSeedPoint[2] + fSeedDirection[2] - SPxyz[2],
                SPxyz[0] - fSeedPoint[0],
                SPxyz[1] - fSeedPoint[1],
                SPxyz[2] - fSeedPoint[2],
                crossprod[0],
                crossprod[1],
                crossprod[2]);

      return pow(pow(crossprod[0], 2) + pow(crossprod[1], 2) + pow(crossprod[2], 2), 0.5) /
             pow(pow(fSeedDirection[0], 2) + pow(fSeedDirection[1], 2) + pow(fSeedDirection[2], 2),
                 0.5);
    }
  }

  //----------------------------------------------------------------------------
  int Seed::GetPointingSign(Seed const& AnotherSeed) const
  {
    double OtherPos[3], OtherErr[3];
    AnotherSeed.GetPoint(OtherPos, OtherErr);
    double DotProd = (OtherPos[0] - fSeedPoint[0]) * fSeedDirection[0] +
                     (OtherPos[1] - fSeedPoint[1]) * fSeedDirection[1] +
                     (OtherPos[2] - fSeedPoint[2]) * fSeedDirection[2];
    return ((DotProd > 0) - (DotProd < 0));
  }

  //----------------------------------------------------------------------------
  void CrossProd(double x1,
                 double x2,
                 double x3,
                 double y1,
                 double y2,
                 double y3,
                 double& out1,
                 double& out2,
                 double& out3)
  {
    out1 = (x2 * y3 - x3 * y2);
    out2 = (x3 * y1 - x1 * y3);
    out3 = (x1 * y2 - x2 * y1);
  }

  //----------------------------------------------------------------------
  std::ostream& operator<<(std::ostream& o, Seed const& a)
  {
    o << "Printing seed contents : " << a.fSeedPoint[0] << " " << a.fSeedPoint[1] << " "
      << a.fSeedPoint[2] << ", " << a.fSeedDirection[0] << " " << a.fSeedDirection[1] << " "
      << a.fSeedDirection[2];

    return o;
  }

  //----------------------------------------------------------------------
  // < operator.
  //
  bool operator<(const Seed& a, const Seed& b)
  {
    if (a.fSeedPoint[2] != b.fSeedPoint[2]) { return a.fSeedPoint[2] < b.fSeedPoint[2]; }
    else if (a.fSeedPoint[1] != b.fSeedPoint[1]) {
      return a.fSeedPoint[1] < b.fSeedPoint[1];
    }

    return a.fSeedPoint[0] < b.fSeedPoint[0];
  }

}
