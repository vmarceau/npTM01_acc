/**
 * @file    paraxialTM01.hpp
 * @brief   Header file that contains functions pertaining to gaussian paraxial TM01 pulses

 * This header file contains the functions pertaining to gaussian paraxial TM01 pulses propagating in free space. This
   includes the closed form expressions of the electromagnetic fields, as well as systems of ODEs describing the motion of
   charged particles in the electromagnetic fields.
 * @author  Vincent Marceau (vincent.marceau.2@ulaval.ca)
 * @since   November 2011
 * @date    November 2011
 * @see C. Varin, "Impulsions d'électrons relativistes ultrarapides à l'aide d'un schéma d'accélération par laser dans le vide,"
   PhD Thesis, Université Laval, 2006.
 * @see P.-L. Fortin, "Dynamique d'un nuage d'électrons soumis à un faisceau TM01 ultra-intense
   et ultrabref," MSc Thesis, Université Laval, 2008.
 */

#ifndef PTM01_HPP_INCLUDED
#define PTM01_HPP_INCLUDED

#include <cmath>
#include <complex>
#include <vector>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>

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
 * @brief Parameter structure for paraxial gaussian TM01 pulses
 */
struct Pparams {
  /** Amplitude parameter [V/m] */
  const double E0;
  /** Rayleigh distance parameter [m] */
  const double zR;
  /** Pulse duration [s] */
  const double T;
  /** Angular frequency [rad/s] */
  const double omega0;
  /** Wave vector [1/m] */
  const double k0;
  /** Phase [rad] */
  const double phi0;
};


/**
 * @brief Computes the electric and magnetic fields of a paraxial TM01 pulse under the form of an analytic signal
 *
 * The radial and longitudinal electric field components are computed from the formulas
 * @f[ E_r = E_0 e^{1/2} \left( \frac{j z_R}{\tilde{q}(z)} \right)^2 \left( \frac{k_0}{z_R} \right)^{1/2} r
  e^{ - jk_0 r^2/2\tilde{q}(z)} e^{- t'^2/T^2} e^{j(\omega_0 t - k_0 z - \phi_0 )} \ ,  @f]
 * and
 * @f[ E_z = E_0 e^{1/2} \left( \frac{j z_R}{\tilde{q}(z)} \right)^2 \left( \frac{k_0}{z_R} \right)^{1/2} \left( -\frac{2j}{k_0}
  \right) \left( 1- \frac{jk_0 r^2}{2\tilde{q}(z)} \right) e^{ - jk_0 r^2/2\tilde{q}(z)} e^{- t'^2/T^2}
   e^{j(\omega_0 t - k_0 z - \phi_0 )} \ ,  \ .@f]
 * where @f$ t' \equiv t - z/c @f$. The azimutal magnetic field component is computed from the formula
   @f[ H_\phi =  \frac{E_r}{\eta_0}  \ . @f]
 * IMPORTANT: The physical electromagnetic fields of the pulse correspond to the real part of the analytic signal.
 * @param E0 Amplitude parameter \f$ E_0 \f$
 * @param zR Rayleigh distance \f$  z_R \f$
 * @param T Duration of the pulse \f$ T \f$
 * @param omega0 Frequency \f$ \omega_0 \f$
 * @param phi0 Phase \f$ \phi_0 \f$
 * @param r Radial coordinate \f$ r \f$
 * @param z Longitudinal coordinate \f$ z \f$
 * @param t Time \f$ t \f$
 * @return Three dimensional vector \f$ (E_r, E_z, H_\phi) \f$
 */
vector<complex<double> > get_pTM01fields(const double E0, const double zR, const double T, const double omega0,
                              const double phi0, const double r, const double z, const double t)
{

  // Various quantities that will be necessary in the computation of the electric field components
  double k0 = omega0/C;
  complex<double> qtilde = z+I*zR;
  complex<double> amplitude = E0*exp(0.5)*pow(I*zR/qtilde,2.0)*sqrt(k0/zR)*exp(-I*k0*pow(r,2.0)/(2.0*qtilde));
  complex<double> pulse_shape = exp(-pow(t-z/C,2.0)/pow(T,2.0));
  complex<double> propagator = exp(I*(omega0*t - k0*z - phi0));
  complex<double> common_factor = amplitude*pulse_shape*propagator;

  // Radial electric field component
  complex<double> Er = common_factor*r;

  // Longitudinal electric field component
  complex<double> Ez = common_factor*(-2.0*I/k0)*(1.0-I*k0*pow(r,2.0)/(2.0*qtilde));

  // Azimutal magnetic field component
  complex<double> Hp = EPS0*Er*C;

  // Return the components of the electric field
  vector<complex<double> > TM01fields;
  TM01fields.push_back(Er);
  TM01fields.push_back(Ez);
  TM01fields.push_back(Hp);
  return TM01fields;

} // end function get_pTM01field


/**
 * @brief System of ODEs describing the motion of a single electron on the optical axis in the electromagnetic fields of paraxial
   gaussian TM01 pulse.
 *
 * The equations of motion of a single electron on the optical axis in the electromagnetic fields of a paraxial gaussian TM01
   pulse are
 * @f[ \frac{dz}{dt} = v_z \ , @f]
 * @f[ \frac{dv_z}{dt} = - \frac{e}{m_e} E_z(z,t) \left( 1- \frac{v_z^2}{c^2} \right)^{3/2} \ , @f]
 * with \f$ z(0) = 0 \f$ , \f$ v_z(0) = 0 \f$ .
 *
 * Since \f$ E_r = H_\phi = 0 \f$ at \f$ r=0 \f$ , the longitudinal electric field is directly computed inside this function from
   the expression
 * @f[ E_z = E_0 e^{1/2} \left( \frac{j z_R}{\tilde{q}(z)} \right)^2 \left( \frac{k_0}{z_R} \right)^{1/2} \left( -\frac{2j}{k_0}
   \right) e^{- t'^2/T^2} e^{j(\omega_0 t - k_0 z - \phi_0 )} \ ,  \ .@f]
 * where @f$ t' \equiv t - z/c @f$ , @f$ \tilde{q}(z) \equiv z +j z_R @f$.
 *
 * @param t Current time \f$ t \f$
 * @param y[] Pointer to the position vector \f$  y=(z,v_z) \f$
 * @param f[] Pointer to a container that will hold \f$ dy/dt \f$
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return GSL_SUCCESS
 */
int dydt_pTM01_onaxis(double t, const double y[], double f[], void * param) {

  // Cast parameters
  Pparams& p = *static_cast<Pparams* >(param);

  // Create multi_array references
  typedef boost::multi_array_ref<const double,1> CSTDBLvecref;
  typedef boost::multi_array_ref<double,1> DBLvecref;
  typedef CSTDBLvecref::index index;
  CSTDBLvecref::extent_gen extents1;
  DBLvecref::extent_gen extents2;
  CSTDBLvecref Y(y,extents1[2]);
  DBLvecref F(f,extents2[2]);

  // Compute the longitudinal electric field
  complex<double> qtilde = Y[0]+I*p.zR;
  double tprime = t - Y[0]/C;
  double Ez = real( -2.0*I*p.E0*exp(0.5)*pow(I*p.zR/qtilde,2.0)*pow(p.zR*p.k0,-0.5)*exp(-pow(tprime/p.T,2.0))*exp(I*(p.omega0*t-p.k0*Y[0]-p.phi0)) );

  // Equations of motion:
  F[0] = Y[1]; // dzdt
  F[1] = -(QE/ME)*Ez*pow(1-pow(Y[1],2.0)/pow(C,2.0),1.5); // dvz/dt

  return GSL_SUCCESS;

} // end function dydt_pTM01_onaxis

#endif // PTM01_HPP_INCLUDED
