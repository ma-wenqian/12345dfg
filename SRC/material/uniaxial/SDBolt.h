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
// $Source: ./OpenSees/SRC/material/uniaxial/SDBolt.h,v $
                                                                        
                                                                        
#ifndef SDBolt_h
#define SDBolt_h

// Written: WXG 
// Created: 
// Revision: 
//
// Description: This file contains the implementation of the joint bolt model   
// that considers contact


#include <UniaxialMaterial.h>
#define MAT_TAG_SDBolt 20240518
// Default values for slops Kse, Kl, Ksp, and the lock threshold 
//#define BOLT_DEFAULT_Kse             50.0
//#define BOLT_DEFAULT_Kl              50.0
//#define BOLT_DEFAULT_Ksp             50.0
#define BOLT_DEFAULT_Threshold       50.0

class SDBolt : public UniaxialMaterial
{
  public:
    SDBolt(int tag, double fsy, double Kse, double Kl, double Ksp, double Threshd = BOLT_DEFAULT_Threshold);
    SDBolt();
    ~SDBolt();

    const char *getClassType(void) const {return "SDBolt";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain(void);              
    double getStress(void);
    double getTangent(void);
    double getInitialTangent(void) {return Kse * Kl / (Kse + Kl);};

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
// AddingSensitivity:BEGIN //////////////////////////////////////////
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////
	//by SAJalali
	virtual double getEnergy() { return Energy; }

 protected:
    
 private:
	/*** Material Properties ***/
    double fy;  // Yield stress
    double Kse;  // Stiffness Kse
    double Kl;   // Stiffness Kl (will go to a large value after the threshold) 
    double Ksp;  // After yielding stiffness 
    double Threshd; //The lock threshold
    
    double Kl_init;// The initial Kl
    double Kll; // Kl aftering locking

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

  double Energy;
  
    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);

    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
