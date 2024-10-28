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
                                                                        
// This file uses support from OpenAI.                                                                        
// Written: Wenqian Ma
// Created: 2024-10
// Revision: 
//
// Description: This file contains the class implementation for 
// Joint01. 
//
// What: "@(#) Joint01.C, revA"

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
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 20) {
    opserr << "Invalid #args,  want: uniaxialMaterial Joint01 tag? Totally 20 parameters are required. K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,\n K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,\n Kl, fsy, Kse, Ksp, Fdp, Fdn " << endln;
    return 0;
  }

  // 这个函数用来读取输入的参数
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args,  want: uniaxialMaterial Joint01 tag? Totally 20 parameters are required. K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,\n K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,\n Kl, fsy, Kse, Ksp, Fdp, Fdn " << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new Joint01(iData[0], dData[0], dData[1], dData[2],dData[3], dData[4], dData[5], dData[6], dData[7], 
                                      dData[8], dData[9], dData[10], dData[11], dData[12], dData[13], dData[14], 
                                      dData[15], dData[16], dData[17], dData[18], dData[19]);
  if (theMaterial == 0) {
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
    Kl(Kl), fsy(fsy), Kse(Kse), Ksp(Ksp), Fdp(Fdp), Fdn(Fdn) {
    trialStrain = 0.0;
    trialStrainRate = 0.0;
    G3p=G2p; G3n=G2n;
    Vsy=fsy / (Kse * Kl / (Kse + Kl));
    V1yp=G2n;
    V2yp=V1yp+fbyp/K2ep;
    V3yp=V1yp+fbyp/K3ep;
    Vdp=V1yp+(Fdp-fbyp)/K1pp;
    V1yn=-G2n;
    V2yn=V1yn-fbyn/K2en;
    V3yn=V1yn-fbyn/K3en;
    Vdn=V1yn-(Fdn-fbyn)/K1pn;
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

}


Joint01::Joint01()
  : UniaxialMaterial(0, MAT_TAG_Joint01), K1ep(0.0), fbyp(0.0), K1pp(0.0), K2ep(0.0), K3ep(0.0), G1p(0.0), G2p(0.0),
    K1en(0.0), fbyn(0.0), K1pn(0.0), K2en(0.0), K3en(0.0), G1n(0.0), G2n(0.0),
    Kl(0.0), fsy(0.0), Kse(0.0), Ksp(0.0), Fdp(0.0), Fdn(0.0) {
    trialStrain = 0.0;
    trialStrainRate = 0.0;
}


Joint01::~Joint01()
{
  // does nothing
}

int Joint01::setTrialStrain(double strain, double strainRate)
{
  trialStrain = strain;
  trialStrainRate = strainRate;

  // 计算每种材料的应力
  trialStress = 0.0;
  trialTangent = 0.0;
  calculate_Bolt(&trialStress, &trialTangent);
  calculate_Glubam1(&trialStress, &trialTangent);
  calculate_Glubam2(&trialStress, &trialTangent);
  calculate_Glubam3(&trialStress, &trialTangent);
  return 0;
}

double Joint01::getStrain(void)
{
  return trialStrain;
}

double Joint01::getStress(void)
{
  return trialStress;
}

double Joint01::getTangent(void)
{
  return trialTangent;
}

double Joint01::getInitialTangent(void)
{
  return Kse * Kl / (Kse + Kl);
}

int Joint01::commitState(void)
{
  return 0;
}

int Joint01::revertToLastCommit(void)
{
  return 0;
}

int Joint01::revertToStart(void)
{
  trialStrain = 0.0;
  trialStrainRate = 0.0;
  return 0;
}

UniaxialMaterial *Joint01::getCopy(void) {
    return new Joint01(this->getTag(), K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,
                       K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,
                       Kl, fsy, Kse, Ksp, Fdp, Fdn);
}

// int 
// Joint01::sendSelf(int cTag, Channel &theChannel)
// {
//   int res = 0;
//   static Vector data(5);
//   data(0) = this->getTag();
//   data(1) = Epos;
//   data(2) = Eneg;
//   data(3) = eta;
//   data(4) = parameterID;
//   res = theChannel.sendVector(this->getDbTag(), cTag, data);
//   if (res < 0) 
//     opserr << "Joint01::sendSelf() - failed to send data" << endln;

//   return res;
// }


// int 
// Joint01::recvSelf(int cTag, Channel &theChannel, 
// 			  FEM_ObjectBroker &theBroker)
// {
//   int res = 0;
//   static Vector data(5);
//   res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
//   if (res < 0) {
//     opserr << "Joint01::recvSelf() - failed to receive data" << endln;
//     Epos = Eneg = 0; 
//     this->setTag(0);      
//   }
//   else {
//     this->setTag(int(data(0)));
//     Epos = data(1);
//     Eneg = data(2);
//     eta  = data(3);
//     parameterID = (int)data(4);
//   }
    
//   return res;
// }

int Joint01::sendSelf(int commitTag, Channel &theChannel) {
    // 发送自身状态
    return 0;
}

int Joint01::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    // 接收自身状态
    return 0;
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

void Joint01::calculate_Bolt(double *trialStress, double *trialTangent)
{
  double v = trialStrain;
  if (v < -Vsy)
  {
    *trialStress += Ksp * (v + Vsy) - fsy;
    *trialTangent += Ksp;
  }
  else if (v < Vsy)
  {
    *trialStress += Kse * Kl / (Kse + Kl) * v;
    *trialTangent += Kse * Kl / (Kse + Kl);
  }
  else
  {
    *trialStress += Ksp * (v - Vsy) + fsy;
    *trialTangent += Ksp;
  }
}

void Joint01::calculate_Glubam1(double *trialStress, double *trialTangent)
{
  double v = trialStrain;

  if (v < Vdn)
  {
    *trialStress += -0.8 * Fdn;
    *trialTangent += 0.0;
  }
  else if (v < V1yn)
  {
    *trialStress += K1pn * (v - V1yn) - fbyn;
    *trialTangent += K1pn;
  }

  else if (v < -G1n)
  {
    *trialStress += K1en * (v + G1n);
    *trialTangent += K1en;
  }
  else if (v < G1p)
  {
    *trialStress += 0.0;
    *trialTangent += 0.0;
  }
  else if (v < V1yp)
  {
    *trialStress += K1ep * (v - V1yp);
    *trialTangent += K1ep;
  }
  else if (v < Vdp)
  {
    *trialStress += K1pp * (v - G2p) + fbyp;
    *trialTangent += K1pp;
  }
  else
  {
    *trialStress += 0.8 * Fdp;
    *trialTangent += 0.0;
  }
}

void Joint01::calculate_Glubam2(double *trialStress, double *trialTangent)
{
  double v = trialStrain;

  if (v < V2yn)
  {
    *trialStress += -fbyn;
    *trialTangent += 0.0;
  }
  else if (v < V1yn)
  {
    *trialStress += K2en * (v - V1yn);
    *trialTangent += K2en;
  }
  else if (v < V1yp)
  {
    *trialStress += 0.0;
    *trialTangent += 0.0;
  }
  else if (v < V2yp)
  {
    *trialStress += K2ep * (v - V1yp);
    *trialTangent += K2ep;
  }
  else
  {
    *trialStress += fbyp;
    *trialTangent += 0.0;
  }
}

void Joint01::calculate_Glubam3(double *trialStress, double *trialTangent)
{
  double v = trialStrain;

  if (v < V3yn)
  {
    *trialStress += -fbyn;
    *trialTangent += 0.0;
  }
  else if (v < V1yn)
  {
    *trialStress += K3en * (v - V1yn);
    *trialTangent += K3en;
  }
  else if (v < V1yp)
  {
    *trialStress += 0.0;
    *trialTangent += 0.0;
  }
  else if (v < V3yp)
  {
    *trialStress += K3ep * (v - V1yp);
    *trialTangent += K3ep;
  }
  else
  {
    *trialStress += fbyp;
    *trialTangent += 0.0;
  }
}