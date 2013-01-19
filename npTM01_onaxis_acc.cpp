/**
 * @file    npTM01_onaxis_acc.cpp
 * @brief   Source file used to simulate the acceleration of a single electron on the optical axis by an ultrashort nonparaxial
            TM01 pulse.

 * Source file used to simulate the acceleration of a single electron on the optical axis by an ultrashort nonparaxial
   TM01 pulse. Results are also compared with their couterparts that use the corresponding paraxial TM01 pulse.
 * @author  Vincent Marceau (vincent.marceau.2@ulaval.ca)
 * @since   November 2011
 * @date    November 2011
 * @see C. Varin, "Impulsions d'électrons relativistes ultrarapides à l'aide d'un schéma d'accélération par laser dans le vide,"
   PhD Thesis, Université Laval, 2006.
 * @see P.-L. Fortin, "Dynamique d'un nuage d'électrons soumis à un faisceau TM01 ultra-intense
   et ultrabref," MSc Thesis, Université Laval, 2008.
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include <boost/multi_array.hpp>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

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
#ifndef ME
/** @brief Electron mass (in kg) */
#define ME GSL_CONST_MKSA_MASS_ELECTRON
#endif // ME
#ifndef QE
/** @brief Electron charge (in C) */
#define QE GSL_CONST_MKSA_ELECTRON_CHARGE
#endif // QE


using namespace std;

/**
 * @brief Main function
 */
int main()
{

  // PARAMETERS

  // Nonparaxial beam
  const double E0 = 3.32096e13;  // Amplitude parameter [V/m]
  const double lambda0 = 800e-9; // Wavelength [m]
  const double k0 = 2.0*PI/lambda0; // Wave vector [1/m]
  const double a = 494.479/k0; // Confocal parameter [m]
  const double omega0 = k0*C; // Angular frequency [rad/s]
  const double s = 277.0; // Spectral width parameter []
  const double phi0 = 0.0; // Pulse phase [rad]
  const double sigma_t = 1.0*s/(omega0*sqrt(2.0*s-1.0)); // Electric energy density temporal standard deviation [S]
  vector<complex<double> > G0 = get_G0_npTM01(phi0,omega0,s); // Amplitudes of the G_\pm^n factors [1/s^n]

  NPparams npparam = {E0,a,s,omega0,G0}; // Parameter container

  // Paraxial beam (additional parameters)
  //const double zR = a; // Rayleigh distance [m]
  const double zR = 6.2832e-5;
  //const double T = sqrt(2.0*s)/omega0; // Pulse length [s]
  const double T = 10e-15;

  Pparams pparam = {E0,zR,T,omega0,k0,phi0}; // Parameter container

  // Integrator parameters
  double t = -100e-15;
  double dt = 1e-17;
  double t_step = 5e-16;
  double t_max = 10e-12;
  const double eps_abs = 1e-20;
  const double eps_rel = 4e-12;

  // Define GSL odeiv parameters
  const gsl_odeiv_step_type * step_type = gsl_odeiv_step_rkf45; // Runge-Kutta-Fehlberg 4-5 stepper
  gsl_odeiv_step * step = gsl_odeiv_step_alloc (step_type,2);
  gsl_odeiv_control * control = gsl_odeiv_control_y_new (eps_abs,eps_rel);
  gsl_odeiv_evolve * evolve = gsl_odeiv_evolve_alloc (2);
  gsl_odeiv_system sys = {dydt_npTM01_onaxis, NULL, 2, &npparam};

  // Open output files
  ofstream out1;
  out1.open("./dat/test/npTM01_test.dat", ios::out);
  out1 << "# t     z     v_z     W" << endl;
  ofstream out2;
  out2.open("./dat/test/pTM01_test.dat", ios::out);
  out2 << "# t     z     v_z     W" << endl;


  // NUMERICAL INTEGRATION

  // Create position and velocity vector
  typedef boost::multi_array<double,1> DBLvec;
  typedef DBLvec::index index;
  DBLvec::extent_gen extents;
  DBLvec y(extents[2]);

  // 1. Acceleration using the nonparaxial pulse

  // Set initial conditions
  y[0] = -lambda0; // z(0)
  y[1] = 0.0; // vz (0)

  // Numerical integration for nonparaxial pulse
  int status = GSL_SUCCESS;
  double t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
    out1 << t << "  " << y[0] << "  " << y[1] << "  " << (1.0-pow(1-pow(y[1],2.0)/pow(C,2.0),-0.5))*ME*pow(C,2.0)/(1.0e6*QE) << endl;
  } // end for


  // 2. Acceleration using the paraxial pulse

  // Change system of equations to be integrated
  sys = {dydt_pTM01_onaxis, NULL, 2, &pparam};

  // Reset initial conditions
  y[0] = -lambda0; // z(0)
  y[1] = 0.0; // vz (0)

  // Numerical integration for paraxial pulse
  t = -100e-15;
  dt = 1e-17;
  status = GSL_SUCCESS;
  t_target = t;
  for (; t_target <= t_max; t_target += t_step ) {
    while (t < t_target) {
      status = gsl_odeiv_evolve_apply (evolve,control,step,&sys,&t,t_target,&dt,y.data());
      if (status != GSL_SUCCESS)
        break;
    } // end while
    if (status != GSL_SUCCESS)
      break;
    out2 << t << "  " << y[0] << "  " << y[1] << "  " << pow(1-pow(y[1],2.0)/pow(C,2.0),-0.5)*ME*pow(C,2.0)/(1.0e6*QE) << endl;
  } // end for


  return 0;
}
