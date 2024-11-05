/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.0 $
// $Date: 2024-10-15 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Joint01.cpp,v $
                                                                        
                                                                      
// Written: Wenqian Ma
// Created: 2024-10
// Revision: 
//
// Description: This file contains the class implementation for 
// Joint01. 
//
// What: "@(#) Joint01.C, revA"


// #define dbl_Epsilon 2.22e-16	 // a very small double 
// already defined this project in other material files

#include <Joint01.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <math.h>
#include <float.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void *
OPS_Joint01(void)
{
  // print out some KUDO's
  opserr << "Joint01 unaxial material - Written by Wenqian Ma, 2024-10-15" << endln;
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int iData[1];
  double dData[20];
  int numData = 1;
  // 这个函数用来读取输入的material tag
  if (OPS_GetIntInput(&numData, iData) != 0)
  {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 20)
  {
    opserr << "Invalid #args,  want: uniaxialMaterial Joint01 tag? Totally 20 parameters are required. K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,\n K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,\n Kl, fsy, Kse, Ksp, Fdp, Fdn " << endln;
    return 0;
  }

  // 这个函数用来读取输入的参数
  if (OPS_GetDoubleInput(&numData, dData) != 0)
  {
    opserr << "Invalid #args,  want: uniaxialMaterial Joint01 tag? Totally 20 parameters are required. K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,\n K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,\n Kl, fsy, Kse, Ksp, Fdp, Fdn " << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new Joint01(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
                            dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14],
                            dData[15], dData[16], dData[17], dData[18], dData[19]);
  if (theMaterial == 0)
  {
    opserr << "WARNING could not create uniaxialMaterial of type Joint01" << endln;
    return 0;
  }
  return theMaterial;
}

Joint01::Joint01(int tag, double K1ep, double fbyp, double K1pp, double K2ep, double K3ep, double G1p, double G2p,
                 double K1en, double fbyn, double K1pn, double K2en, double K3en, double G1n, double G2n,
                 double Kl, double fsy, double Kse, double Ksp, double Fdp, double Fdn)
    : UniaxialMaterial(tag, MAT_TAG_Joint01), K1ep(K1ep), fbyp(fbyp), K1pp(K1pp), K2ep(K2ep), K3ep(K3ep), G1p(G1p), G2p(G2p),
      K1en(K1en), fbyn(fbyn), K1pn(K1pn), K2en(K2en), K3en(K3en), G1n(G1n), G2n(G2n),
      Kl(Kl), fsy(fsy), Kse(Kse), Ksp(Ksp), Fdp(Fdp), Fdn(Fdn)
{
  this->revertToStart();
}

Joint01::Joint01()
    : UniaxialMaterial(0, MAT_TAG_Joint01), K1ep(0.0), fbyp(0.0), K1pp(0.0), K2ep(0.0), K3ep(0.0), G1p(0.0), G2p(0.0),
      K1en(0.0), fbyn(0.0), K1pn(0.0), K2en(0.0), K3en(0.0), G1n(0.0), G2n(0.0),
      Kl(0.0), fsy(0.0), Kse(0.0), Ksp(0.0), Fdp(0.0), Fdn(0.0)
{
  this->revertToStart();
}

Joint01::~Joint01()
{
  // does nothing
}

int Joint01::setTrialStrain(double strain, double strainRate)
{

  trialStrainRate = strainRate;

  // 计算每种材料的应力
  trialStress = 0.0;
  trialTangent = 0.0;
  // Reset history variables to last converged state
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  Tstress_Bolt = Cstress_Bolt;
  Tstress_Glubam1p = Cstress_Glubam1p;
  Tstress_Glubam1n = Cstress_Glubam1n;
  Tstress_Glubam2p = Cstress_Glubam2p;
  Tstress_Glubam2n = Cstress_Glubam2n;
  Tstress_Glubam3p = Cstress_Glubam3p;
  Tstress_Glubam3n = Cstress_Glubam3n;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;
  // opserr <<"strain: "<< strain << "\tCstrain: " << Cstrain << "\tdStrain: " << dStrain << endln;

  if (fabs(dStrain) > DBL_EPSILON)
  {
    // Set trial strain
    Tstrain = strain;

    calculate_Bolt(dStrain, &Tstress_Bolt, &Cstress_Bolt, &trialTangent);
    calculate_Glubam(&Tstrain, &Tstress_Glubam1p, &trialTangent,
                     &minElasticYieldStrain1p, &maxElasticYieldStrain1p, fbyp, K1ep, G1p, K1pp, &Vdp);
    calculate_Glubam(&Tstrain, &Tstress_Glubam1n, &trialTangent,
                     &minElasticYieldStrain1n, &maxElasticYieldStrain1n, -fbyn, K1en, -G1n, K1pn, &Vdn);
    calculate_Glubam(&Tstrain, &Tstress_Glubam2p, &trialTangent,
                     &minElasticYieldStrain2p, &maxElasticYieldStrain2p, fbyp, K2ep, G2p);
    calculate_Glubam(&Tstrain, &Tstress_Glubam2n, &trialTangent,
                     &minElasticYieldStrain2n, &maxElasticYieldStrain2n, -fbyn, K2en, -G2n);
    calculate_Glubam(&Tstrain, &Tstress_Glubam3p, &trialTangent,
                     &minElasticYieldStrain3p, &maxElasticYieldStrain3p, fbyp, K3ep, G3p);
    calculate_Glubam(&Tstrain, &Tstress_Glubam3n, &trialTangent,
                       &minElasticYieldStrain3n, &maxElasticYieldStrain3n, -fbyn, K3en, -G3n);

    Tstress = Tstress_Bolt + Tstress_Glubam1p+ Tstress_Glubam1n + Tstress_Glubam2p + Tstress_Glubam2n + Tstress_Glubam3p + Tstress_Glubam3n;
  }
  return 0;
}

double Joint01::getStrain(void)
{
  return Tstrain;
}

double Joint01::getStress(void)
{
  return Tstress;
}

double Joint01::getTangent(void)
{
  return Ttangent;
}

double Joint01::getInitialTangent(void)
{
  return Kse * Kl / (Kse + Kl) + 0 + 0 + 0 + 0;
}

int Joint01::commitState(void)
{
  Cstrain = Tstrain;
  Cstress = Tstress;
  Cstress_Bolt = Tstress_Bolt;
  Cstress_Glubam1p = Tstress_Glubam1p;
  Cstress_Glubam1n = Tstress_Glubam1n;
  Cstress_Glubam2p = Tstress_Glubam2p;
  Cstress_Glubam2n = Tstress_Glubam2n;
  Cstress_Glubam3p = Tstress_Glubam3p;
  Cstress_Glubam3n = Tstress_Glubam3n;
  Ctangent = Ttangent;
  return 0;
}

int Joint01::revertToLastCommit(void)
{
  // Reset trial state variables to last committed state
  Tstrain = Cstrain;
  Tstress = Cstress;
  Tstress_Bolt = Cstress_Bolt;
  Tstress_Glubam1p = Cstress_Glubam1p;
  Tstress_Glubam1n = Cstress_Glubam1n;
  Tstress_Glubam2p = Cstress_Glubam2p;
  Tstress_Glubam2n = Cstress_Glubam2n;
  Tstress_Glubam3p = Cstress_Glubam3n;
  Tstress_Glubam3n = Cstress_Glubam3n;
  Ttangent = Ctangent;
  return 0;
}

int Joint01::revertToStart(void)
{
  // caculate the variables
  trialStrain = 0.0;
  trialStrainRate = 0.0;
  K = Kse * Kl / (Kse + Kl);
  G3p = G2p;
  G3n = G2n;
  Vsy = fsy / (Kse * Kl / (Kse + Kl));
  V1yp = G1p+fbyp/K1ep;
  V2yp = G2p + fbyp / K2ep;
  V3yp = G3p + fbyp / K3ep;
  Vdp = V1yp + (Fdp - fbyp) / K1pp;
  V1yn = -G1n - fbyn / K1en;
  V2yn = -G2n - fbyn / K2en;
  V3yn = -G3n - fbyn / K3en;
  Vdn = V1yn - (Fdn - fbyn) / K1pn;

  // Sets all history and state variables to initial values
  // History variables
  Cstrain = 0.0;
  Cstress = 0.0;
  Cstress_Bolt = 0.0;
  Cstress_Glubam1p = 0.0;
  Cstress_Glubam1n = 0.0;
  Cstress_Glubam2p = 0.0;
  Cstress_Glubam2n = 0.0;
  Cstress_Glubam3p = 0.0;
  Cstress_Glubam3n = 0.0;
  Ctangent = Kse * Kl / (Kse + Kl);

  Tstrain = 0.0;
  Tstress = 0.0;
  Tstress_Bolt = 0.0;
  Tstress_Glubam1p = 0.0;
  Tstress_Glubam1n = 0.0;
  Tstress_Glubam2p = 0.0;
  Tstress_Glubam2n = 0.0;
  Tstress_Glubam3p = 0.0;
  Tstress_Glubam3n = 0.0;
  Ttangent = Kse * Kl / (Kse + Kl);

  // set initial values for the glumbam variables
  minElasticYieldStrain1p = G1p;
  maxElasticYieldStrain1p = V1yp;
  minElasticYieldStrain1n = -G1n; // G1n is positive, here should be a negative number  
  maxElasticYieldStrain1n = V1yn;

  minElasticYieldStrain2p = G2p;
  maxElasticYieldStrain2p = V2yp;
  minElasticYieldStrain2n = -G2n; // G2n is positive, here should be a negative number  
  maxElasticYieldStrain2n = V2yn;

  minElasticYieldStrain3p = G3p;
  maxElasticYieldStrain3p = V3yp;
  minElasticYieldStrain3n = -G3n; // G3n is positive, here should be a negative number  
  maxElasticYieldStrain3n = V3yn;

  // 我自己额外增加的 为了看参数有没有顺利传递
  opserr << "Joint01 tag: " << this->getTag() << endln;
  opserr << "Input parameters: " << endln;
  opserr << " K1ep: " << K1ep << " fbyp: " << fbyp << " K1pp: " << K1pp << " K2ep: " << K2ep << " K3ep: " << K3ep << " G1p: " << G1p << " G2p: " << G2p << endln;
  opserr << " K1en: " << K1en << " fbyn: " << fbyn << " K1pn: " << K1pn << " K2en: " << K2en << " K3en: " << K3en << " G1n: " << G1n << " G2n: " << G2n << endln;
  opserr << " Kl: " << Kl << " fsy: " << fsy << " Kse: " << Kse << " Ksp: " << Ksp << " Fdp: " << Fdp << " Fdn: " << Fdn << endln;
  opserr << "Caculated parameters: " << endln;
  opserr << " Vsy: " << Vsy << endln;
  opserr << " V1yp: " << V1yp << " V2yp: " << V2yp << " V3yp: " << V3yp << " Vdp: " << Vdp << endln;
  opserr << " V1yn: " << V1yn << " V2yn: " << V2yn << " V3yn: " << V3yn << " Vdn: " << Vdn << endln;
  return 0;
}

UniaxialMaterial *Joint01::getCopy(void) {
    return new Joint01(this->getTag(), K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,
                       K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,
                       Kl, fsy, Kse, Ksp, Fdp, Fdn);
}

// There are two methods provided which are required is the user uses to use the database or parallel
// procesing features of the OpenSees applications. If neither are to be used, the developer need simply
// return a negative value in both methods. The idea is that the material must pack up it's information
// using Vector and ID objects and send it off to a Channel object. On the flip side, the receiving
// blank element must receive the same Vector and ID data, unpack it and set the variables.

int Joint01::sendSelf(int commitTag, Channel &theChannel) {
    // 发送自身状态
    return -1;
}

int Joint01::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    // 接收自身状态
    return -1;
}

void Joint01::Print(OPS_Stream &s, int flag)
{
  s << "Joint01 tag: " << this->getTag() << endln;
  s << "Input parameters: " << endln;
  s << " K1ep: " << K1ep << " fbyp: " << fbyp << " K1pp: " << K1pp << " K2ep: " << K2ep << " K3ep: " << K3ep << " G1p: " << G1p << " G2p: " << G2p << endln;
  s << " K1en: " << K1en << " fbyn: " << fbyn << " K1pn: " << K1pn << " K2en: " << K2en << " K3en: " << K3en << " G1n: " << G1n << " G2n: " << G2n << endln;
  s << " Kl: " << Kl << " fsy: " << fsy << " Kse: " << Kse << " Ksp: " << Ksp << " Fdp: " << Fdp << " Fdn: " << Fdn << endln;
  s << "Caculated parameters: " << endln;
  s << " Vsy:" << Vsy << endln;
  s << " V1yp:" << V1yp << " V2yp:" << V2yp << " V3yp:" << V3yp << " Vdp:" << Vdp << endln;
  s << " V1yn:" << V1yn << " V2yn:" << V2yn << " V3yn:" << V3yn << " Vdn:" << Vdn << endln;
}

void Joint01::calculate_Bolt(double dStrain, double *Tstress_Bolt, double *Cstress_Bolt, double *trialTangent)
{
  if (fabs(Tstrain) > Vsy)
  {
    K = Kse;
  }
  // hardenring type : kinematic hardening
  // double Y1 = Ksp * Tstrain + fsy * (1.0 - Ksp / (Kse * Kl / (Kse + Kl)));
  // double Y2 = *Cstress_Bolt + Kse * Kl / (Kse + Kl) * dStrain;
  // double Y3 = Ksp * Tstrain - fsy * (1.0 - Ksp / (Kse * Kl / (Kse + Kl)));

  double Y1 = Ksp * Tstrain + fsy - Ksp * Vsy;
  double Y2 = *Cstress_Bolt + K * dStrain;
  double Y3 = Ksp * Tstrain - (fsy - Ksp * Vsy);


  *Tstress_Bolt = std::max(Y3, std::min(Y1, Y2));

  if (*Tstress_Bolt == Y2)
  {
    *trialTangent += Ksp;
  }
  else
  {

    *trialTangent += Kse * Kl / (Kse + Kl);
  }
}


// reference: /OpenSees/SRC/material/uniaxial/EPPGapMaterial.cpp line 140 - 172, 202 - 236
void Joint01::calculate_Glubam(double *Tstrain, double *Tstress_Glubam, double *trialTangent,
                               double *minElasticYieldStrain, double *maxElasticYieldStrain, double fy, double E, double gap, double eta_E, double *Vdrop)
{
  // hardenring type : isotropic hardening

  // set the trial strain
  trialStrain = *Tstrain;

  if (Vdrop == nullptr)
  {
    Vdrop = new double(DBL_MAX);
  }


  // determine trial stress and tangent
  if (fy >= 0)
  {
    if (trialStrain > *Vdrop)
    {
      *Tstress_Glubam = 0.8 * fy;
      *trialTangent += 0;
      *maxElasticYieldStrain = DBL_MAX;
      *minElasticYieldStrain = DBL_MAX; // 脆性破坏 brittle failure.
      // *minElasticYieldStrain = trialStrain - *Tstress_Glubam / E; // 塑性破坏 plastic failure.
      *Vdrop = trialStrain;
      return; // 提前退出函数
    }
    else if (trialStrain > *maxElasticYieldStrain)
    {
      *Tstress_Glubam = fy + (trialStrain - gap - fy / E) * eta_E;
      *trialTangent += eta_E;
      *maxElasticYieldStrain = trialStrain;
      *minElasticYieldStrain = trialStrain - *Tstress_Glubam / E;
    }
    else if (trialStrain < *minElasticYieldStrain)
    {
      *Tstress_Glubam = 0;
      *trialTangent += 0;
    }
    else
    {
      *Tstress_Glubam = E * (trialStrain - *minElasticYieldStrain);
      *trialTangent += E;
    }

  }
  else
  {
     if (trialStrain < -fabs(*Vdrop))
    {
      *Tstress_Glubam = 0.8 * fy;
      *trialTangent += 0;
      *maxElasticYieldStrain = -DBL_MAX;
      *minElasticYieldStrain = -DBL_MAX;
      *Vdrop = trialStrain;
      return; // 提前退出函数
    }
    else if (trialStrain < *maxElasticYieldStrain)
    {
      *Tstress_Glubam = fy + (trialStrain - gap - fy / E) * eta_E;
      *trialTangent += eta_E;
      *maxElasticYieldStrain = trialStrain;
      *minElasticYieldStrain = trialStrain - *Tstress_Glubam / E;
    }
    else if (trialStrain > *minElasticYieldStrain)
    {
      *Tstress_Glubam = 0;
      *trialTangent += 0;
    }
    else
    {
      *Tstress_Glubam = E * (trialStrain - *minElasticYieldStrain);
      *trialTangent += E;
    }
  }

}
