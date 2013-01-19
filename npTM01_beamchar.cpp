/**
 * @file    npTM01_beamchar.cpp
 * @brief   Source file used to characterize an ultrashort nonparaxial TM01 beam

 * This source file is used to characterize the various properties of an ultrashort nonparaxial TM01 beam. We also compare the
   nonparaxial beam with its paraxial counterpart.
 * @author  Vincent Marceau (vincent.marceau.2@ulaval.ca)
 * @since   November 2011
 * @date    November 2011
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>

#include "./nonparaxialTM01.hpp"
#include "./paraxialTM01.hpp"


#ifndef I
/** @brief Imaginary unit */
#define I complex<double>(0,1.0)
#endif // I
#ifndef PI
/** @brief Pi */
#define PI M_PI
#endif // PI
#ifndef C
/** @brief Velocity of light in free space (in m/s) */
#define C GSL_CONST_MKSA_SPEED_OF_LIGHT
#endif // C
#ifndef EPS0
/** @brief Permittivity of free space (in F/m) */
#define EPS0 GSL_CONST_MKSA_VACUUM_PERMITTIVITY
#endif // EPS0


using namespace std;

/**
 * @brief Main function
 */
int main() {

  // BEAM PARAMETERS

  // Nonparaxial beam
  const double E0 = 1e10;  // Amplitude parameter [V/m]
  const double lambda0 = 800e-9; // Wavelength [m]
  const double k0 = 2.0*PI/lambda0; // Wave vector [1/m]
  const double a = 1.0/k0; // Confocal parameter [m]
  const double omega0 = k0*C; // Angular frequency [rad/s]
  const double s = 1.0; // Spectral width parameter [] **MUST BE A POSITIVE INTEGER**
  const double T = sqrt(2.0*s)/omega0; // Pulse length [s]
  const double phi0 = 0.0; // Pulse phase [rad]
  const double sigma_t = 1.0*s/(omega0*sqrt(2.0*s-1.0)); // Electric energy density temporal standard deviation [s]
	const double w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a,2.0))-1)/k0; // Beam waist size [m]
	const double zR = k0*pow(w0,2.0)/2; // Rayleigh distance [m]
  vector<complex<double> > G0 = get_G0_npTM01(phi0,omega0,s); // Amplitudes of the G_\pm^n factors [1/s^n]

  // NONPARAXIAL PULSE CHARACTERIZATION

  // Density of electric energy versus radial coordinate r and time t for an ultrashort nonparaxial TM01 pulse
  /*ofstream out1; // Open output file
  out1.open("./dat/beamchar/npTM01_we_vs_rt.dat", ios::out);
  out1 << "# Electric field energy density at beam waist versus radial coordinate and time for an ultrashort nonparaxial TM01
           pulse" << endl;
  out1 << "# Beam parameters:" << endl;
  out1 << "# E0 = " << E0 << " , ";
  out1 << "omega0 = " << omega0 << " , ";
  out1 << "lambda0 = " << lambda0 << " , ";
  out1 << "ka = " << k0*a << endl;
  out1 << "# s = " << s << " , ";
  out1 << "sigma_t = " << sigma_t << " , ";
  out1 << "phi0 = " << phi0 << endl;
  out1 << "#" << endl;
  out1 << "# r/lambda0     t/sigma_t     w_e" << endl;

  for (double r=-5.0*lambda0; r<=5.0*lambda0; r+=lambda0/50.0) {
    for (double t=-5.0*sigma_t; t<=5.0*sigma_t; t+=sigma_t/50.0) {
      vector<complex<double> > TM01fields = get_npTM01fields(E0,a,s,omega0,G0,abs(r),0.0,t); // EM fields
      double w_e = EPS0*(pow(abs(TM01fields[0]),2.0) + pow(abs(TM01fields[1]),2.0))/2.0;
      out1 << r/lambda0 << "  " << t/sigma_t << "  " << w_e << endl;
    } // end for
    out1 << endl;
  } // end for*/


  // Instantaneous power along the optical axis versus time for an ultrashort nonparaxial TM01 pulse
  /* ofstream out2; // Open output file
  out2.open("./dat/beamchar/npTM01_p_vs_t.dat", ios::out);
  out2 << "# Instantaneous power in the z direction at beam waist versus time for an ultrashort nonparaxial TM01 pulse" << endl;
  out2 << "# Beam parameters:" << endl;
  out2 << "# E0 = " << E0 << " , ";
  out2 << "omega0 = " << omega0 << " , ";
  out2 << "lambda0 = " << lambda0 << " , ";
  out2 << "ka = " << k0*a << endl;
  out2 << "# s = " << s << " , ";
  out2 << "sigma_t = " << sigma_t << " , ";
  out2 << "phi0 = " << phi0 << endl;
  out2 << "#" << endl;
  out2 << "# t     P_z" << endl;

  double dr = lambda0/200.0;
  for (double t=-5.0*sigma_t; t<=5.0*sigma_t; t+=sigma_t/50.0) {
    double power = 0.0;
    for (double r=0; r<=30.0*lambda0; r+=dr) {
      vector<complex<double> > TM01fields = get_npTM01fields(E0,a,s,omega0,G0,r,0.0,t); // EM fields
      double Sz = real(TM01fields[0])*real(TM01fields[2]) ; // Poynting vector z component (Er*Hphi)
      power += Sz*(2.0*PI*r*dr + PI*dr*dr);
    } // end for
    out2 << t << "  " << power << endl;
  } // end for*/


  // Peak power along the optical axis versus confocal parameter for an ultrashort nonparaxial TM01 pulse
  ofstream out3; // Open output file
  out3.open("./dat/2012.01.26.beamchar/npTM01_p_vs_a.dat", ios::out);
  out3 << "# Peak power in the z direction at beam waist versus confocal parameter for an ultrashort nonparaxial TM01 pulse";
  out3 << endl;
  out3 << "# Beam parameters:" << endl;
  out3 << "# E0 = " << E0 << " , ";
  out3 << "omega0 = " << omega0 << " , ";
  out3 << "lambda0 = " << lambda0 << " , ";
  out3 << "# s = " << s << " , ";
  out3 << "sigma_t = " << sigma_t << " , ";
  out3 << "phi0 = " << phi0 << endl;
  out3 << "#" << endl;
  out3 << "# a     P_z (manual)     P_z (gsl)" << endl;

  double dr = lambda0/500.0;
  for (double avar=1.0/k0; avar<=100.0/k0; avar+=1.0/k0) {

    // Manual integration
    double power_man = 0.0;
    double w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*avar,2.0))-1)/k0;
    for (double r=0; r<=10.0*w0; r+=dr) {
      vector<complex<double> > TM01fields = get_npTM01fields(E0,avar,s,omega0,G0,r,0.0,0.0); // EM fields
      double Sz = real(TM01fields[0])*real(TM01fields[2]) ; // Poynting vector z component (Er*Hphi)
      power_man += Sz*(2.0*PI*r*dr + PI*dr*dr);
    } // end for

    // Integration with GSL
    NPparams npparam = {E0,avar,s,omega0,G0};
    double power_gsl = get_Ppeak_npTM01 (10.0*w0,&npparam);
    vector<complex<double> > TM01fields = get_npTM01fields(E0,avar,s,omega0,G0,0.0,0.0,0.0);

    out3 << k0*avar << "  " << power_man << "  " << power_gsl << "  " << abs(TM01fields[1]) << endl;
  } // end for


  // Peak power along the optical axis versus pulse duration parameter for an ultrashort nonparaxial TM01 pulse
  /*ofstream out4; // Open output file
  out4.open("./dat/beamchar/npTM01_p_vs_s.dat", ios::out);
  out4 << "# Peak power in the z direction at beam waist versus pulse duration for an ultrashort nonparaxial TM01 pulse" << endl;
  out4 << "# Beam parameters:" << endl;
  out4 << "# E0 = " << E0 << " , ";
  out4 << "omega0 = " << omega0 << " , ";
  out4 << "lambda0 = " << lambda0 << " , ";
  out4 << "ka = " << k0*a << endl;
  out4 << "# sigma_t = " << sigma_t << " , ";
  out4 << "phi0 = " << phi0 << endl;
  out4 << "#" << endl;
  out4 << "# s     P_z (manual)     P_z (gsl)" << endl;

  double dr = lambda0/1000.0;
  double w0 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a,2.0))-1)/k0;

  for (double svar=1.0; svar<=1000.0; svar+=2.0) {

    G0 = get_G0_npTM01(phi0,omega0,svar);

    // Manual integration
    double power_man = 0.0;
    for (double r=0; r<=10.0*w0; r+=dr) {
      vector<complex<double> > TM01fields = get_npTM01fields(E0,a,svar,omega0,G0,r,0.0,0.0); // EM fields
      double Sz = real(TM01fields[0])*real(TM01fields[2]) ; // Poynting vector z component (Er*Hphi)
      power_man += Sz*(2.0*PI*r*dr + PI*dr*dr);
    } // end for

    // Integration with GSL
    NPparams npparam = {E0,a,svar,omega0,G0};
    double power_gsl = get_Ppeak_npTM01 (10.0*w0,&npparam);

    out4 << svar << "  " << power_man << "  " << power_gsl << endl;

  } // end for */


  // Peak power along the optical axis versus E0 for an ultrashort nonparaxial TM01 pulse
  /*ofstream out9; // Open output file
  out9.open("./dat/beamchar/npTM01_p_vs_E0.dat", ios::out);
  out9 << "# Peak power in the z direction at beam waist amplitude parameter E0 for an ultrashort nonparaxial TM01 pulse" << endl;
  out9 << "# Beam parameters:" << endl;
  out9 << "# omega0 = " << omega0 << " , ";
  out9 << "lambda0 = " << lambda0 << " , ";
  out9 << "phi0 = " << phi0 << endl;
  out9 << "#" << endl;
  out9 << "# E0     a=1,s=1     a=50,s=1     a=1,s=100     a=50,s=100" << endl;

  double a1 = 1.0/k0;
  double a2 = 50.0/k0;
  double s1 = 1.0;
  double s2 = 100.0;

  double w01 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a1,2.0))-1)/k0;
  double w02 = sqrt(2.0)*sqrt(sqrt(1.0+pow(k0*a2,2.0))-1)/k0;
  vector<complex<double> > G01 = get_G0_npTM01(phi0,omega0,s1);
  vector<complex<double> > G02 = get_G0_npTM01(phi0,omega0,s2);

  for (double E0var=1.0; E0var<=1000.0; E0var+=1.0) {

    // Integration with GSL
    NPparams npparam1 = {E0var,a1,s1,omega0,G01};
    NPparams npparam2 = {E0var,a2,s1,omega0,G01};
    NPparams npparam3 = {E0var,a1,s2,omega0,G02};
    NPparams npparam4 = {E0var,a2,s2,omega0,G02};

    double power_gsl1 = get_Ppeak_npTM01 (10.0*w01,&npparam1);
    double power_gsl2 = get_Ppeak_npTM01 (10.0*w02,&npparam2);
    double power_gsl3 = get_Ppeak_npTM01 (10.0*w01,&npparam3);
    double power_gsl4 = get_Ppeak_npTM01 (10.0*w02,&npparam4);

    out9 << E0var << "  " << power_gsl1 << "  " << power_gsl2 << "  " << power_gsl3 << "  " << power_gsl4 << endl;

  } // end for*/


  // COMPARISON BETWEEN NONPARAXIAL AND PARAXIAL PULSES

  // Peak power along the optical axis versus confocal parameter for the paraxial versus the nonparaxial TM01 pulse
  /*ofstream out5; // Open output file
  out5.open("./dat/beamchar/pnpTM01_p_vs_a.dat", ios::out);
  out5 << "# Peak power in the z direction at beam waist versus confocal parameter for a paraxial versus a nonparaxial TM01 pulse"
       << endl;
  out5 << "# Beam parameters of the nonparaxial pulse:" << endl;
  out5 << "# E0 = " << E0 << " , ";
  out5 << "omega0 = " << omega0 << " , ";
  out5 << "lambda0 = " << lambda0 << " , ";
  out5 << "# s = " << s << " , ";
  out5 << "sigma_t = " << sigma_t << " , ";
  out5 << "phi0 = " << phi0 << endl;
  out5 << "# Additional beam parameters for the paraxial pulse:" << endl;
  out5 << "# T = " << T << endl;
  out5 << "#" << endl;
  out5 << "# ka     P_z (npar)     P_z (par)" << endl;

  double dr = lambda0/200.0;
  for (double avar=100.0/k0; avar<=1000.0/k0; avar+=10.0/k0) {
    double powerNP = 0.0;
    double powerP = 0.0;
    for (double r=0; r<=45.0*lambda0; r+=dr) {
      vector<complex<double> > npTM01fields = get_npTM01fields(E0,avar,s,omega0,G0,r,0.0,0.0); // Nonparaxial pulse
      vector<complex<double> > pTM01fields = get_pTM01fields(E0,avar,T,omega0,phi0,r,0.0,0.0); // Paraxial pulse
      double SzNP = real(npTM01fields[0])*real(npTM01fields[2]); // Poynting vector z component (nonparaxial)
      double SzP = real(pTM01fields[0])*real(pTM01fields[2]) ; // Poynting vector z component (paraxial)
      powerNP += SzNP*(2.0*PI*r*dr + PI*dr*dr);
      powerP += SzP*(2.0*PI*r*dr + PI*dr*dr);
    } // end for
    out5 << k0*avar << "  " << powerNP << "  " << powerP << endl;
  } // end for/*


  // On-axis longitudinal electric field versus z at t=0 for the paraxial versus the nonparaxial TM01 pulse
  /*ofstream out6; // Open output file
  out6.open("./dat/beamchar/pnpTM01_Ez_vs_z.dat", ios::out);
  out6 << "# On-axis longitudinal electric field versus z at t=0 for the paraxial versus the nonparaxial TM01 pulse" << endl;
  out6 << "# Beam parameters of the nonparaxial pulse:" << endl;
  out6 << "# E0 = " << E0 << " , ";
  out6 << "omega0 = " << omega0 << " , ";
  out6 << "lambda0 = " << lambda0 << " , ";
  out6 << "ka = " << k0*a << endl;
  out6 << "# s = " << s << " , ";
  out6 << "sigma_t = " << sigma_t << " , ";
  out6 << "phi0 = " << phi0 << endl;
  out6 << "# Additional beam parameters for the paraxial pulse:" << endl;
  out6 << "# kzR = " << k0*zR << " , ";
  out6 << "T = " << T << endl;
  out6 << "#" << endl;
  out6 << "# z/lambda0     E_z (npar)     E_z (par)" << endl;

  for (double z=-5.0*lambda0; z<=5.0*lambda0; z+=lambda0/50.0) {
    vector<complex<double> > npTM01fields = get_npTM01fields(E0,a,s,omega0,G0,0.0,z,0.0); // Nonparaxial pulse
    vector<complex<double> > pTM01fields = get_pTM01fields(E0,zR,T,omega0,phi0,0.0,z,0.0); // Paraxial pulse
    out6 << z/lambda0 << "  " << real(npTM01fields[1]) << "  " << real(pTM01fields[1]) << endl;
  } // end for*/


  // Ezmax/Ermax versus confocal parameter for the paraxial versus the nonparaxial TM01 pulse
  /*ofstream out7; // Open output file
  out7.open("./dat/beamchar/pnpTM01_Ez_vs_a.dat", ios::out);
  out7 << "# Ezmax/Ermax versus confocal parameter for a paraxial versus a nonparaxial TM01 pulse"
      << endl;
  out7 << "# Beam parameters of the nonparaxial pulse:" << endl;
  out7 << "# E0 = " << E0 << " , ";
  out7 << "omega0 = " << omega0 << " , ";
  out7 << "lambda0 = " << lambda0 << " , ";
  out7 << "# s = " << s << " , ";
  out7 << "sigma_t = " << sigma_t << " , ";
  out7 << "phi0 = " << phi0 << endl;
  out7 << "# Additional beam parameters for the paraxial pulse:" << endl;
  out7 << "# T = " << T << endl;
  out7 << "#" << endl;
  out7 << "# ka     Ezmax/Ermax (npar)     Ezmax/Ermax (par)" << endl;

  for (double avar=0.1/k0; avar<=25.0/k0; avar+=0.05/k0) {
    // Peak value of the longitudinal electric field is at r=0, z=0, t=0;
    vector<complex<double> > npTM01fields = get_npTM01fields(E0,avar,s,omega0,G0,0.0,0.0,0.0); // Nonparaxial pulse
    double npEzmax = abs(npTM01fields[1]);
    vector<complex<double> > pTM01fields = get_pTM01fields(E0,avar,T,omega0,phi0,0.0,0.0,0.0); // Paraxial pulse
    double pEzmax = abs(pTM01fields[1]);
    // Find peak value of the radial electric field
    double npErmax = 0.0;
    double pErmax = 0.0;
    double dr = lambda0/500.0;
    for (double r=0; r<=30.0*lambda0; r+=dr) {
      vector<complex<double> > npTM01fields = get_npTM01fields(E0,avar,s,omega0,G0,r,0.0,0.0); // Nonparaxial pulse
      vector<complex<double> > pTM01fields = get_pTM01fields(E0,avar,T,omega0,phi0,r,0.0,0.0); // Paraxial pulse
      if (abs(npTM01fields[0]) > npErmax) npErmax = abs(npTM01fields[0]);
      if (abs(pTM01fields[0]) > pErmax) pErmax = abs(pTM01fields[0]);
    } // end for
    out7 << k0*avar << "  " << npEzmax/npErmax << "  " << pEzmax/pErmax << endl;
  } // end for*/


  // Electric field at beam waist versus r at t=0 for the paraxial versus the nonparaxial TM01 pulse
  /*ofstream out8; // Open output file
  out8.open("./dat/beamchar/pnpTM01_E_vs_r.dat", ios::out);
  out8 << "# Electric field at beam waist versus r at t=0 for the paraxial versus the nonparaxial TM01 pulse" << endl;
  out8 << "# Beam parameters of the nonparaxial pulse:" << endl;
  out8 << "# E0 = " << E0 << " , ";
  out8 << "omega0 = " << omega0 << " , ";
  out8 << "lambda0 = " << lambda0 << " , ";
  out8 << "ka = " << k0*a << endl;
  out8 << "# s = " << s << " , ";
  out8 << "sigma_t = " << sigma_t << " , ";
  out8 << "phi0 = " << phi0 << endl;
  out8 << "# Additional beam parameters for the paraxial pulse:" << endl;
  out8 << "# kzR = " << k0*zR << " , ";
  out8 << "T = " << T << endl;
  out8 << "#" << endl;
  out8 << "# r/lambda0     abs(E_r) (npar)     abs(E_z) (npar)     abs(E_r) (par)     abs(E_z) (par)" << endl;

  double dr = lambda0/50.0;
  for (double r=0.0; r<=20000.0*lambda0; r+=dr) {
    vector<complex<double> > npTM01fields = get_npTM01fields(E0,a,s,omega0,G0,r,0.0,0.0); // Nonparaxial pulse
    vector<complex<double> > pTM01fields = get_pTM01fields(E0,zR,T,omega0,phi0,r,0.0,0.0); // Paraxial pulse
    //out8 << r/lambda0 << "  " << abs(npTM01fields[0]) << "  " << abs(npTM01fields[1])
    //                  << "  " << abs(pTM01fields[0]) << "  " << abs(pTM01fields[1]) << endl;
    out8 << r/lambda0 << "  " << r*real(npTM01fields[0])*real(npTM01fields[2])
                      << "  " << r*real(pTM01fields[0])*real(pTM01fields[2]) << endl;
  } // end for */


  // Theoretical energy gain limit for paraxial and nonparaxial TM01 pulses
	/*double ElimP_sum = 0.0;
	double ElimNP_sum = 0.0;
	double ElimNP_gsl = 0.0;

  // Integration by direct summation
	double dz=zR/100.0;
	double b=5000.0*zR;
  for (double z=0; z<=b; z+=dz) {
    vector<complex<double> > npTM01fields = get_npTM01fields(E0,a,s,omega0,G0,0.0,z,z/C); // Nonparaxial pulse
    vector<complex<double> > pTM01fields = get_pTM01fields(E0,zR,T,omega0,phi0,0.0,z,z/C); // Paraxial pulse
		ElimNP_sum+=real(npTM01fields[1])*dz;
		ElimP_sum+=real(pTM01fields[1])*dz;
  } // end for

  // Integration with gsl
  NPparams npparam = {E0,a,s,omega0,G0};
  void * param = &npparam;
  double error;
  int key = 6; // 61 points Gauss-Kronrod rule
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(20000);
  gsl_function F;
  F.function = &get_Ezforward_npTM01;
  F.params = param;
  // Perform numerical integration using GSL QAG adaptive integration for smooth functions
  gsl_integration_qag(&F,0,b,0,1e-6,20000,key,w,&ElimNP_gsl,&error);
  // Free allocated workspace
  gsl_integration_workspace_free(w);

  cout << "Paraxial beam (summation): Elim = " << 1.0e-6*ElimP_sum << " MeV" << endl;
  cout << "Nonparaxial beam (summation): Elim = " << 1.0e-6*ElimNP_sum << " MeV" << endl;
  cout << "Nonparaxial beam (GSL integration): Elim = " << 1.0e-6*ElimNP_gsl << " MeV, error = " << 1.0e-6*error << endl;*/

  return 0;
}
