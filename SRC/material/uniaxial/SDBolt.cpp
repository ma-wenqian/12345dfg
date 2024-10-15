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
                                                                        
// $Revision: 0.10 $
// $Date: 2024-05-18 $
// $Source: ./OpenSees/SRC/material/uniaxial/SDBolt.cpp,v $
                                                                        
// Written: WXG 
// Created: 
// Revision: 
//
// Description: This file contains the implementation of the joint bolt model   
// that considers contact


#include <SDBolt.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <math.h>
#include <float.h>


#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_SDBolt()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[5];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial SDBolt tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 4 && numData != 5) {
    opserr << "Invalid #args, want: uniaxialMaterial SDBolt " << iData[0] << " fy? Kse? Kl? Ksp? <Lock Threshold?>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial SDBolt " << iData[0] << " fy? Kse? Kl? Ksp? <Lock Threshold?>" << endln;
    return 0;
  }

  if (numData == 4) {
    dData[4] = BOLT_DEFAULT_Threshold;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SDBolt(iData[0], dData[0], dData[1],
						   dData[2], dData[3], dData[4]);

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type SDBolt Material" << endln;
    return 0;
  }

  return theMaterial;
}



SDBolt::SDBolt
(int tag, double FY, double KSE, double KL, double KSP, double THRESHD):
   UniaxialMaterial(tag, 0),
   fy(FY), Kse(KSE), Kl(KL), Ksp(KSP), Threshd(THRESHD),
   Energy(0.0), parameterID(0), SHVs(0)
{  
   //Set Kl inital value and Kll value
   Kl_init = Kl;
   Kll = Kse * 9999; //set Kll to a large value

   // Sets all history and state variables to initial values
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Kse * Kl / (Kse + Kl);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Kse * Kl / (Kse + Kl);
}

SDBolt::SDBolt():UniaxialMaterial(0,0),
		   fy(0.0), Kse(0.0), Kl(0.0), Ksp(0.0), Threshd(0.0),
		   Energy(0.0), parameterID(0), SHVs(0)
{

}

SDBolt::~SDBolt()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int SDBolt::setTrialStrain (double strain, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   return 0;
}

int SDBolt::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
   // Reset history variables to last converged state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   // Determine change in strain from last converged state
   double dStrain = strain - Cstrain;

   if (fabs(dStrain) > DBL_EPSILON) {
     // Set trial strain
     Tstrain = strain;

     // Calculate the trial state given the trial strain
     determineTrialState (dStrain);

   }

   stress = Tstress;
   tangent = Ttangent;

   return 0;
}

void SDBolt::determineTrialState (double dStrain)
{
      double fyOneMinusB = fy * (1.0 - Ksp * (Kse + Kl) / Kse / Kl); // where b in steel01 is Ksp * (Kse + Kl) / Kse / Kl

      double Esh = Ksp;
	  double epsy = fy * (Kse + Kl) / Kse / Kl ;
      
      double c1 = Esh*Tstrain;
      
      double c2 = TshiftN*fyOneMinusB;

      double c3 = TshiftP*fyOneMinusB;

      double c = Cstress + Kse * Kl / (Kse + Kl) *dStrain;

      /**********************************************************
         removal of the following lines due to problems with
	 optimization may be required (e.g. on gnucc compiler
         with optimization turned on & -ffloat-store option not
         used) .. replace them with line that follows but which 
         now requires 2 function calls to achieve same result !!
      ************************************************************/

      double c1c3 = c1 + c3;

      if (c1c3 < c)
	Tstress = c1c3;
      else
	Tstress = c;

      double c1c2 = c1-c2;

      if (c1c2 > Tstress)
	Tstress = c1c2;

      /* ***********************************************************
      and replace them with:

      Tstress = fmax((c1-c2), fmin((c1+c3),c));
      **************************************************************/

      if (fabs(Tstress-c) < DBL_EPSILON)
	  Ttangent = Kse * Kl / (Kse + Kl);
      else
	Ttangent = Esh;

      //
      // Determine if a load reversal has occurred due to the trial strain
      //

      // Determine initial loading condition
      if (Tloading == 0 && dStrain != 0.0) {
	  if (dStrain > 0.0)
	    Tloading = 1;
	  else
	    Tloading = -1;
      }

      // Transition from loading to unloading, i.e. positive strain increment
      // to negative strain increment
      if (Tloading == 1 && dStrain < 0.0) {
	  Tloading = -1;
	  if (Cstrain > TmaxStrain)
	    TmaxStrain = Cstrain;
	  //TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
	  TshiftN = 1 + pow((TmaxStrain - TminStrain) / (2.0 * epsy), 0.8);
	  }

      // Transition from unloading to loading, i.e. negative strain increment
      // to positive strain increment
      if (Tloading == -1 && dStrain > 0.0) {
	  Tloading = 1;
	  if (Cstrain < TminStrain)
	    TminStrain = Cstrain;
	  //TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
	  TshiftP = 1 + pow((TmaxStrain - TminStrain) / (2.0 * epsy), 0.8);
      }
}

void SDBolt::detectLoadReversal (double dStrain)
{
   // Determine initial loading condition
   if (Tloading == 0 && dStrain != 0.0)
   {
      if (dStrain > 0.0)
         Tloading = 1;
      else
         Tloading = -1;
   }

   double epsy = fy/ (Kse * Kl / (Kse + Kl));

   // Transition from loading to unloading, i.e. positive strain increment
   // to negative strain increment
   if (Tloading == 1 && dStrain < 0.0)
   {
      Tloading = -1;
      if (Cstrain > TmaxStrain)
         TmaxStrain = Cstrain;
      //TshiftN = 1 + a1*pow((TmaxStrain-TminStrain)/(2.0*a2*epsy),0.8);
	  TshiftN = 1 + pow((TmaxStrain - TminStrain) / (2.0 * epsy), 0.8);
   }

   // Transition from unloading to loading, i.e. negative strain increment
   // to positive strain increment
   if (Tloading == -1 && dStrain > 0.0)
   {
      Tloading = 1;
      if (Cstrain < TminStrain)
         TminStrain = Cstrain;
      //TshiftP = 1 + a3*pow((TmaxStrain-TminStrain)/(2.0*a4*epsy),0.8);
	  TshiftP = 1 + pow((TmaxStrain - TminStrain) / (2.0 * epsy), 0.8);
   }
}

double SDBolt::getStrain ()
{
   return Tstrain;
}

double SDBolt::getStress ()
{
   return Tstress;
}

double SDBolt::getTangent ()
{
   return Ttangent;
}

int SDBolt::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CmaxStrain = TmaxStrain;
   CshiftP = TshiftP;
   CshiftN = TshiftN;
   Cloading = Tloading;

   // State variables
   //by SAJalali
   Energy += 0.5*(Tstress + Cstress)*(Tstrain - Cstrain);

   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   fabs(Cstrain) > Threshd ? Kl = Kll : Kl = Kl_init;

   return 0;
}

int SDBolt::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;

   // Reset trial state variables to last committed state
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int SDBolt::revertToStart ()
{
   // History variables
   CminStrain = 0.0;
   CmaxStrain = 0.0;
   CshiftP = 1.0;
   CshiftN = 1.0;
   Cloading = 0;

   TminStrain = 0.0;
   TmaxStrain = 0.0;
   TshiftP = 1.0;
   TshiftN = 1.0;
   Tloading = 0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Kse * Kl / (Kse + Kl);

   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Kse * Kl / (Kse + Kl);

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* SDBolt::getCopy ()
{
	SDBolt* theCopy = new SDBolt(this->getTag(), fy, Kse, Kl, Ksp, Threshd);

   // Converged history variables
   theCopy->CminStrain = CminStrain;
   theCopy->CmaxStrain = CmaxStrain;
   theCopy->CshiftP = CshiftP;
   theCopy->CshiftN = CshiftN;
   theCopy->Cloading = Cloading;

   // Trial history variables
   theCopy->TminStrain = TminStrain;
   theCopy->TmaxStrain = TmaxStrain;
   theCopy->TshiftP = TshiftP;
   theCopy->TshiftN = TshiftN;
   theCopy->Tloading = Tloading;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   // Trial state variables
   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   theCopy->parameterID = parameterID;
   if (SHVs != 0)
     theCopy->SHVs = new Matrix(*SHVs);
   
   return theCopy;
}

int SDBolt::sendSelf (int commitTag, Channel& theChannel)
{
   static Vector data(16);
   data(0) = this->getTag();

   // Material properties
   data(1) = fy;
   data(2) = Kse;
   data(3) = Kl;
   data(4) = Ksp;
   data(5) = Threshd;

   // History variables from last converged state
   data(6) = CminStrain;
   data(7) = CmaxStrain;
   data(8) = CshiftP;
   data(9) = CshiftN;
   data(10) = Cloading;

   // State variables from last converged state
   data(11) = Cstrain;
   data(12) = Cstress;
   data(13) = Ctangent;
   data(14) = parameterID;
   data(15) = -1;
   if (SHVs != 0)
     data(15) = SHVs->noCols();
   
   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   int dbTag = this->getDbTag();
   
   if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
     opserr << "SDBolt::sendSelf() - failed to send data" << endln;
     return -1;
   }

   if (SHVs != 0) {
     if (theChannel.sendMatrix(dbTag, commitTag, *SHVs) < 0) {
       opserr << "SDBolt::sendSelf() - failed to send SHVs matrix" << endln;
       return -2;
     }
   }
   
   return 0;
}

int SDBolt::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  
   static Vector data(16);
  
   if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
     opserr << "SDBolt::recvSelf() - failed to receive data" << endln;
      this->setTag(0);
      return -1;
   }

   this->setTag(int(data(0)));

   // Material properties
   fy = data(1);
   Kse = data(2);
   Kl = data(3);
   Ksp = data(4);
   Threshd = data(5);
   
   // History variables from last converged state
   CminStrain = data(6);
   CmaxStrain = data(7);
   CshiftP = data(8);
   CshiftN = data(9);
   Cloading = int(data(10));
   
   // Copy converged history values into trial values since data is only
   // sent (received) after convergence
   TminStrain = CminStrain;
   TmaxStrain = CmaxStrain;
   TshiftP = CshiftP;
   TshiftN = CshiftN;
   Tloading = Cloading;
   
   // State variables from last converged state
   Cstrain = data(11);
   Cstress = data(12);
   Ctangent = data(13);      
   
   // Copy converged state values into trial values
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;
   
   parameterID = (int)data(14);

   // Receive sensitivity history variables (SHVs)
   int noCols = (int)data(15);   
   if (noCols > 0) {
     if (SHVs != 0)
       delete SHVs;
     SHVs = new Matrix(2, noCols);
     if (SHVs == 0) {
       opserr << "SDBolt::recvSelf() - failed to allocate SHVs matrix" << endln;
       return -2;
     }
     
     if (theChannel.recvMatrix(dbTag, commitTag, *SHVs) < 0) {
       opserr << "SDBolt::recvSelf() - failed to receive SHVs matrix" << endln;
       return -3;
     }
   }      
   
   return 0;
}

void SDBolt::Print (OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {    
    s << "SDBolt tag: " << this->getTag() << endln;
    s << "  fy: " << fy << " ";
    s << "  Kse: " << Kse << " ";
    s << "  Kl: " << Kl << " ";
    s << "  Ksp: " << Ksp << " ";
    s << "  Threshold: " << Threshd << " ";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"SDBolt\", ";
	s << "\"fy\": " << fy << ", ";
	s << "\"Kse\": " << Kse << ", ";
	s << "\"Kl\": " << Kl << ", ";
    s << "\"Ksp\": " << Ksp << ", ";
    s << "\"Threshold\": " << Threshd << ", ";
  }
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
SDBolt::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(fy);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Kse") == 0) {
    param.setValue(Kse);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Kl") == 0) {
    param.setValue(Kl);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Ksp") == 0) {
    param.setValue(Ksp);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"Threshd") == 0) {
    param.setValue(Threshd);
    return param.addObject(5, this);
  }

  return -1;
}



int
SDBolt::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->fy = info.theDouble;
		break;
	case 2:
		this->Kse = info.theDouble;
		break;
	case 3:
		this->Kl = info.theDouble;
		break;
	case 4:
		this->Ksp = info.theDouble;
		break;
	case 5:
		this->Threshd= info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = Kse * Kl / (Kse + Kl);          // Initial stiffness

	return 0;
}




int
SDBolt::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}



double
SDBolt::getStressSensitivity(int gradIndex, bool conditional)
{
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + Kse * Kl / (Kse + Kl) *dStrain;
	double fyOneMinusB = fy * (1.0 - Ksp * (Kse + Kl) / Kse / Kl);
	double Esh = Ksp;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity 
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity * (Ksp * (Kse + Kl) / Kse / Kl)*Tstrain
				 + Kse * Kl / (Kse + Kl) *bSensitivity*Tstrain
				 + TshiftP*(fySensitivity*(1- Ksp * (Kse + Kl) / Kse / Kl)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 - Kse * Kl / (Kse + Kl) *CstrainSensitivity;
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity * (Ksp * (Kse + Kl) / Kse / Kl)*Tstrain
			     + Kse * Kl / (Kse + Kl) *bSensitivity*Tstrain
				 - TshiftN*(fySensitivity*(1- Ksp * (Kse + Kl) / Kse / Kl)-fy*bSensitivity);
	}

	return gradient;
}




double
SDBolt::getInitialTangentSensitivity(int gradIndex)
{
	// For now, assume that this is only called for initial stiffness 
	if (parameterID == 2) {
		return 1.0; 
	}
	else {
		return 0.0;
	}
}


int
SDBolt::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize unconditaional stress sensitivity
	double gradient = 0.0;


	// Pick up sensitivity history variables
	double CstrainSensitivity = 0.0;
	double CstressSensitivity	 = 0.0;
	if (SHVs != 0) {
		CstrainSensitivity = (*SHVs)(0,gradIndex);
		CstressSensitivity = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fySensitivity = 0.0;
	double E0Sensitivity = 0.0;
	double bSensitivity = 0.0;
	if (parameterID == 1) {
		fySensitivity = 1.0;
	}
	else if (parameterID == 2) {
		E0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		bSensitivity = 1.0;
	}


	// Compute min and max stress
	double Tstress;
	double dStrain = Tstrain-Cstrain;
	double sigmaElastic = Cstress + Kse * Kl / (Kse + Kl) *dStrain;
	double fyOneMinusB = fy * (1.0 - Ksp * (Kse + Kl) / Kse / Kl);
	double Esh = Ksp;
	double c1 = Esh*Tstrain;
	double c2 = TshiftN*fyOneMinusB;
	double c3 = TshiftP*fyOneMinusB;
	double sigmaMax = c1+c3;
	double sigmaMin = c1-c2;


	// Evaluate stress sensitivity ('gradient')
	if ( (sigmaMax < sigmaElastic) && (fabs(sigmaMax-sigmaElastic)>1e-5) ) {
		Tstress = sigmaMax;
		gradient = E0Sensitivity * (Ksp * (Kse + Kl) / Kse / Kl)*Tstrain
				 + Kse * Kl / (Kse + Kl) *bSensitivity*Tstrain
				 + Ksp*TstrainSensitivity
				 + TshiftP*(fySensitivity*(1- Ksp * (Kse + Kl) / Kse / Kl)-fy*bSensitivity);
	}
	else {
		Tstress = sigmaElastic;
		gradient = CstressSensitivity 
			     + E0Sensitivity*(Tstrain-Cstrain)
				 + Kse * Kl / (Kse + Kl) *(TstrainSensitivity-CstrainSensitivity);
	}
	if (sigmaMin > Tstress) {
		gradient = E0Sensitivity * (Ksp * (Kse + Kl) / Kse / Kl)*Tstrain
			     + Kse * Kl / (Kse + Kl) *bSensitivity*Tstrain
			     + Ksp *TstrainSensitivity
				 - TshiftN*(fySensitivity*(1- Ksp * (Kse + Kl) / Kse / Kl)-fy*bSensitivity);
	}


	// Commit history variables
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////

