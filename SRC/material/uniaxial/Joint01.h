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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Joint01.h,v $

                                                                        
                                                                        
#ifndef Joint01_h
#define Joint01_h


// Written: Wenqian Ma
// Created: 2024-10
// Revision: 

//
// Description: This file contains the class definition for 
// Joint01. Joint01 provides the abstraction
// of an viscoelastic uniaxial material,
// Totally 20 parameters are required.
// K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p,
// K1en, fbyn, K1pn, K2en, K3en, G1n, G2n,
// Kl, fsy, Kse, Ksp, Fdp, Fdn


#include <UniaxialMaterial.h>

#define MAT_TAG_Joint01 20241014

class Joint01 : public UniaxialMaterial
{
  public:
    Joint01(int tag, double K1ep, double fbyp, double K1pp, double K2ep, double K3ep, double G1p, double G2p, 
            double K1en, double fbyn, double K1pn, double K2en, double K3en, double G1n, double G2n, 
            double Kl, double fsy, double Kse, double Ksp, double Fdp, double Fdn);
    
    Joint01();
    ~Joint01();

    const char *getClassType(void) const {return "Joint01";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);        

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, 
        FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);

    
  private:
    /*** Material Properties ***/
    // 20 input parameters
    double K1ep, fbyp, K1pp, K2ep, K3ep, G1p, G2p;
    double K1en, fbyn, K1pn, K2en, K3en, G1n, G2n;
    double Kl, fsy, Kse, Ksp, Fdp, Fdn;
    double trialStrain;
    double trialStrainRate;
    double trialTangent;
    double trialStress;
    // cālculated values
    double G3p, G3n;
    double Vsy;
    double V1yp, V2yp, V3yp, Vdp;
    double V1yn, V2yn, V3yn, Vdn;

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Cstress_Bolt, Cstress_Glubam1p, Cstress_Glubam1n, Cstress_Glubam2p, Cstress_Glubam2n, Cstress_Glubam3p, Cstress_Glubam3n;
    double Ctangent;   

    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Tstress_Bolt, Tstress_Glubam1p, Tstress_Glubam1n, Tstress_Glubam2p, Tstress_Glubam2n,Tstress_Glubam3p, Tstress_Glubam3n;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience 

    /*** Glubam1 ***/
    double minElasticYieldStrain1p, maxElasticYieldStrain1p;
    double minElasticYieldStrain1n, maxElasticYieldStrain1n;
    
    /*** Glubam2 ***/
    double minElasticYieldStrain2p, maxElasticYieldStrain2p;
    double minElasticYieldStrain2n, maxElasticYieldStrain2n;

    /*** Glubam3 ***/
    double minElasticYieldStrain3p, maxElasticYieldStrain3p;
    double minElasticYieldStrain3n, maxElasticYieldStrain3n;

    // 应力计算方法
    void calculate_Bolt(double dStrain, double *Tstress_Bolt, double * Cstress_Bolt, double* trialTangent);
    void calculate_Glubam(double *Tstrain, double *Tstress_Glubam, double *trialTangent,
                          double *minElasticYieldStrain, double *maxElasticYieldStrain, double fy, double E, double gap,
                          double eta_E = 0, double *Vdrop = nullptr);
  
};


#endif

