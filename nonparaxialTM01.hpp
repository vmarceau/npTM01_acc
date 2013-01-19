/**
 * @file    nonparaxialTM01.hpp
 * @brief   Header file that contains functions pertaining to ultrashort nonparaxial TM01 pulses

 * This header file contains the functions pertaining to ultrashort nonparaxial TM01 pulses propagating in free space. This
   includes the closed form expressions of the electromagnetic fields, as well as systems of ODEs describing the motion of
   charged particles in the electromagnetic fields.
 * @author  Vincent Marceau (vincent.marceau.2@ulaval.ca)
 * @since   November 2011
 * @date    February 2012
 * @see A. April, "Impulsions laser ultrabrèves et fortement focalisées dans le vide," PhD Thesis, Université Laval, 2012.
 * @see C. Varin, "Impulsions d'électrons relativistes ultrarapides à l'aide d'un schéma d'accélération par laser dans le vide,"
   PhD Thesis, Université Laval, 2006.
 * @see P.-L. Fortin, "Dynamique d'un nuage d'électrons soumis à un faisceau TM01 ultra-intense
   et ultrabref," MSc Thesis, Université Laval, 2008.
 */

#ifndef NPTM01_HPP_INCLUDED
#define NPTM01_HPP_INCLUDED

#include <cmath>
#include <complex>
#include <vector>

#include <boost/multi_array.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>


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
 * @brief Parameter structure for ultrashort nonparaxial pulses
 */
struct NPparams {
  /** Amplitude parameter [V/m] */
  const double E0;
  /** Confocal parameter [m] */
  const double a;
  /** Spectral width parameter [] **MUST BE A POSITIVE INTEGER** */
  const double s;
  /** Angular frequency [rad/s] */
  const double omega0;
  /** Amplitudes of the G_\pm^n factors [1/s^n] */
  vector<complex<double> > G0;
};


/**
 * @brief Computes the constant factor that appears in the \f$ G_\pm^{(n)} \f$ factors in the expression of the electromagnetic
   fields of an ultrashort nonparaxial TM01 pulse
 *
 * The \f$ G_\pm^{(n)} \f$ factors can be expressed as
 * @f[  e^{-j\phi_0} \frac{\Gamma(s+n+1)}{\Gamma(s+1)} \left( \frac{j \omega_0}{s}\right)^n \left[
   \left(1- \frac{j \omega_0 \tilde{t}_+}{s}\right)^{-(s+n+1)} \pm
   \left(1- \frac{j \omega_0 \tilde{t}_-}{s} \right)^{-(s+n+1)}\right] \ . @f]
 * Here we compute
 * @f[ G_0^{(n)} \equiv e^{-j\phi_0} \frac{\Gamma(s+n+1)}{\Gamma(s+1)} \left( \frac{j \omega_0}{s}\right)^n \ .@f]
 * @param phi0 Phase \f$ \phi_0 \f$
 * @param omega0 Frequency of maximum amplitude \f$ \omega_0 \f$
 * @param s Spectral width parameter \f$ s \f$
 * @return Three dimensional vector \f$ (G_0^{(0)}, G_0^{(1)}, G_0^{(2)}) \f$
 */
vector<complex<double> > get_G0_npTM01(const double phi0, const double omega0, const double s) {

  complex<double> G00 = exp(-I*phi0);
  complex<double> G01 = G00*(s+1.0)*(I*omega0/s);
  complex<double> G02 = G01*(s+2.0)*(I*omega0/s);

  vector<complex<double> > G0;
  G0.push_back(G00);
  G0.push_back(G01);
  G0.push_back(G02);

  return G0;

} // end function get_G0_npTM01


/**
 * @brief Computes the electric and magnetic fields of an ultrashort nonparaxial TM01 pulse under the form of an analytic signal
 *
 * The radial and longitudinal electric field components are computed from the formulas
 * @f[ E_r = -3 E_0 e^{1/2} \left( \frac{a}{k_0} \right)^{3/2} \frac{\cos \tilde{\theta} \sin \tilde{\theta}}{\tilde{R}} \left(
    \frac{G_-^{(0)}}{\tilde{R}^2} - \frac{G_+^{(1)}}{c\tilde{R}} +  \frac{G_-^{(2)}}{3c^2} \right)  @f]
 * and
 * @f[ E_z =  -2 E_0 e^{1/2} \left( \frac{a}{k_0} \right)^{3/2} \frac{1}{\tilde{R}} \left[
    \frac{1}{2} \left( 3\cos^2 \tilde{\theta} - 1 \right) \left(
    \frac{G_-^{(0)}}{\tilde{R}^2} - \frac{G_+^{(1)}}{c\tilde{R}} +  \frac{G_-^{(2)}}{3c^2} \right)  -
    \frac{G_-^{(2)}}{3c^2}\right]  \ .@f]
 * The azimutal magnetic field component is computed from the formula
   @f[ H_\phi =  -E_0 e^{1/2} \left( \frac{a}{k_0} \right)^{3/2} \frac{\sin \tilde{\theta}}{\eta_0 \tilde{R}} \left(
    \frac{G_-^{(1)}}{c\tilde{R}} - \frac{G_+^{(2)}}{c^2}\right)  \ . @f]
 * IMPORTANT: The physical electromagnetic fields of the pulse correspond to the real part of the analytic signal.
 *
 * @param psi0 Amplitude parameter \f$ \Psi_0 \f$
 * @param a Confocal parameter \f$  a \f$
 * @param s Spectral width parameter \f$ s \f$
 * @param omega0 Frequency of maximum amplitude \f$ \omega_0 \f$
 * @param G0 Coefficients of the \f$ G_\pm^{(n)} \f$ functions
 * @param r Radial coordinate \f$ r \f$
 * @param z Longitudinal coordinate \f$ z \f$
 * @param t Time \f$ t \f$
 * @return Three dimensional vector \f$ (E_r, E_z, H_\phi) \f$
 */
vector<complex<double> > get_npTM01fields(const double E0, const double a, const double s, const double omega0,
                              vector<complex<double> >& G0, const double r, const double z, const double t)
{

  // Various quantities that will be necessary in the computation of the electric field components
  double psi0 = -exp(0.5)*pow(a*C/omega0,1.5)*E0;
  complex<double> ztilde = z + I*a;
  complex<double> Rtilde = sqrt(pow(r,2.0) + pow(ztilde,2.0));
  complex<double> Rtildesquare = pow(r,2.0) + pow(ztilde,2.0);
  complex<double> sintheta = r/Rtilde;
  complex<double> costheta = ztilde/Rtilde;
  complex<double> tplus = t + Rtilde/C + I*a/C;
  complex<double> tminus = t - Rtilde/C + I*a/C;
  double fplusr = abs(1.0-I*omega0*tplus/s);
  double fplusphi = arg(1.0-I*omega0*tplus/s);
  double fminusr = abs(1.0-I*omega0*tminus/s);
  double fminusphi = arg(1.0-I*omega0*tminus/s);
  complex<double> Gminus0 = G0[0]*( pow(fplusr,-(s+1.0))*exp(-I*fplusphi*(s+1.0)) - pow(fminusr,-(s+1.0))*exp(-I*fminusphi*(s+1.0)) );
  complex<double> Gplus1 = G0[1]*( pow(fplusr,-(s+2.0))*exp(-I*fplusphi*(s+2.0)) + pow(fminusr,-(s+2.0))*exp(-I*fminusphi*(s+2.0)) );
  complex<double> Gminus1 = G0[1]*( pow(fplusr,-(s+2.0))*exp(-I*fplusphi*(s+2.0)) - pow(fminusr,-(s+2.0))*exp(-I*fminusphi*(s+2.0)) );
  complex<double> Gplus2 = G0[2]*( pow(fplusr,-(s+3.0))*exp(-I*fplusphi*(s+3.0)) + pow(fminusr,-(s+3.0))*exp(-I*fminusphi*(s+3.0)) );
  complex<double> Gminus2 = G0[2]*( pow(fplusr,-(s+3.0))*exp(-I*fplusphi*(s+3.0)) - pow(fminusr,-(s+3.0))*exp(-I*fminusphi*(s+3.0)) );
  complex<double> Gfactor = Gminus0/Rtildesquare - Gplus1/(C*Rtilde) + Gminus2/(3.0*pow(C,2.0));

  // Radial electric field component
  complex<double> Er = 3.0*psi0*costheta*sintheta*Gfactor/Rtilde;

  // Longitudinal electric field component
  complex<double> Ez = (2.0*psi0/Rtilde)*( 0.5*(3.0*pow(costheta,2.0)-1.0)*Gfactor - Gminus2/(3*pow(C,2.0))  );

  // Azimutal magnetic field component
  complex<double> Hp = psi0*sintheta*EPS0*( Gminus1/Rtilde - Gplus2/C )/Rtilde;

  // Return the components of the electric field
  vector<complex<double> > TM01fields;
  TM01fields.push_back(Er);
  TM01fields.push_back(Ez);
  TM01fields.push_back(Hp);
  return TM01fields;

} // end function get_npTM01fields


/**
 * @brief Computes the z-component of the Poynting vector times \f$ 2\pi r \f$ at beam waist for an ultrashort nonparaxial TM01 pulse
 *
 * The z-component of the Poynting vector at beam waist is given by
 * @f[ S_z(r,z=0,t=0) = E_r(r,z=0,t=0) \times H_\phi (r,z=0,t=0) \ . @f]
 * The z-component is afterwards multiplied by \f$ 2\pi r \f$ for later integration.
 *
 * @param r Radial distance \f$ r \f$
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return Poynting vector \f$ S_z(r,z=0,t=0) \f$ times \\f$ 2\pi r \f$ at radial position \f$ r \f$
 */
 double get_Szwaist_npTM01(double r, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields =  get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,r,0.0,0.0);

  // Compute Sz
  double Sz = 2.0*PI*r*real(npTM01fields[0])*real(npTM01fields[2]);

  return Sz;

} // end function get_Szwaist_npTM01


/**
 * @brief Computes the z-component of the electric field on the optical axis at \f$ t=z/c\f$ for an ultrashort nonparaxial TM01
 * pulse
 * @see Function get_npTM01fields
 *
 * @param z Distance from beam waist on the optical axis
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return Longitudinal electric field \f$ E_z(r=0,z,t=z/c) \f$
 */
 double get_Ezforward_npTM01(double z, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields =  get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,z,z/C);

  // Compute Ez
  double Ez = real(npTM01fields[1]);

  return Ez;

} // end function get_Ezforward_npTM01


/**
 * @brief Computes the z-component of the electric field on the optical axis at \f$ t=-z/c\f$ for an ultrashort nonparaxial TM01
 * pulse
 * @see Function get_npTM01fields
 *
 * @param z Distance from beam waist on the optical axis
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return Longitudinal electric field \f$ E_z(r=0,z,t=-z/c) \f$
 */
 double get_Ezbackward_npTM01(double z, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields =  get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,z,-z/C);

  // Compute Ez
  double Ez = real(npTM01fields[1]);

  return Ez;

} // end function get_Ezbackward_npTM01


/**
 * @brief Computes the peak power of an ultrashort nonparaxial TM01 pulse
 *
 * The peak power is obtained by integrating the z-component of the Poynting vector at beam waist over the x-y space,
 * @f[ P_\textrm{peak} = 2 \pi \int_0^\infty S_z(r,0,0) r dr \ . @f]
 * The numerical integration is performed using GSL.
 *
 * @param b Upper bound for radial integration
 * @param param Void pointer to the parameters.
 * @return Peak power of the pulse
 */
double get_Ppeak_npTM01(double b, void * param) {

  // Get workspace
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);

  // Define the function to be integrated and the result containers
  double Ppeak, error;
  gsl_function F;
  F.function = &get_Szwaist_npTM01;
  F.params = param;

  // Perform numerical integration using GSL QAGS adaptive integration with known singular points
  gsl_integration_qags(&F,0,b,0,1e-6,10000,w,&Ppeak,&error);

  // Free allocated workspace
  gsl_integration_workspace_free(w);

  // Return peak power of the pulse
  return Ppeak;

} // end function get_Ppeak_npTM01


/**
 * @brief Computes the imaginary part of the longitudinal electric field at pulse peak with \f$ \phi0=0 \f$
 * @see Function get_npTM01fields
 *
 * @param z Distance from beam waist on the optical axis
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return Imaginary part of the longitudinal electric field \f$ \textrm{Im}\{ E_z(r=0,z,t=z/c) \} \f$
 */
double get_ImagEz_npTM01(double z, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields = get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,z,z/C);

  // Compute ImagEz
  double ImagEz = imag(npTM01fields[1]);

  return ImagEz;

} // end function get_ImagEz_npTM01


/**
 * @brief Computes the distance from beam waist at which the Gouy phase shift has varied by an amount \f$ \pi/2 \f$
 * @see Function get_ImagEz_npTM01
 *
 * The distance is computed by solving \f$ \textrm{Im}\{ E_z(r=0,z,t=z/c) \} = 0 \f$ for \f$ \phi0=0 \f$.
 * The root finding procedure is carried out using GSL
 *
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return Distance from beam waist at which the Gouy phase shift has varied by an amount \f$ \pi/2 \f$
 */
double get_zG_npTM01(void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Procedure suggested in the GSL reference manual
  int status;
  int iter = 0;
  int max_iter = 200;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double zG = 0;
  double z_lo = 0.5*p.a, z_hi = 10.0*p.a;
  gsl_function F;
  F.function = &get_ImagEz_npTM01;
  F.params = param;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, z_lo, z_hi);

  do {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    zG = gsl_root_fsolver_root (s);
    z_lo = gsl_root_fsolver_x_lower (s);
    z_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (z_lo, z_hi, 0, 0.000001);
  } while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);

  if (iter == max_iter)
    cout << "WARNING: The root finding algorithm for zG has escaped before reaching convergence..." << endl;

  return zG;

} // end function get_zG_npTM01


/**
 * @brief Computes \f$ |E_z(t)|^2 \f$ at beam waist for an ultrashort and nonparaxial TM01 pulsed beam
 *
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return \f$ |E_z(t)|^2 \f$ at beam waist
 */
double get_t0Ez_npTM01(double t, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields = get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,0.0,t);

  // Compute ImagEz
  double t0Ez = pow(real(npTM01fields[1]),2.0);

  return t0Ez;

} // end function get_t0Ez_npTM01


/**
 * @brief Computes \f$ t |E_z(t)|^2 \f$ at beam waist for an ultrashort and nonparaxial TM01 pulsed beam
 *
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return \f$ t |E_z(t)|^2 \f$ at beam waist
 */
double get_t1Ez_npTM01(double t, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields = get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,0.0,t);

  // Compute ImagEz
  double t1Ez = t*pow(real(npTM01fields[1]),2.0);

  return t1Ez;

} // end function get_t1Ez_npTM01


/**
 * @brief Computes \f$ t^2 |E_z(t)|^2 \f$ at beam waist for an ultrashort and nonparaxial TM01 pulsed beam
 *
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return \f$ t^2 |E_z(t)|^2 \f$ at beam waist
 */
double get_t2Ez_npTM01(double t, void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Get electromagnetic fields
  vector<complex<double> > npTM01fields = get_npTM01fields(p.E0, p.a, p.s, p.omega0,p.G0,0.0,0.0,t);

  // Compute ImagEz
  double t2Ez = pow(t,2.0)*pow(real(npTM01fields[1]),2.0);

  return t2Ez;

} // end function get_t2Ez_npTM01


/**
 * @brief System of ODEs describing the motion of a single electron on the optical axis in the electromagnetic fields of an
   ultrashort nonparaxial TM01 pulse.
 *
 * The equations of motion of a single electron on the optical axis in the electromagnetic fields of an ultrashort nonparaxial TM01
   pulse are
 * @f[ \frac{dz}{dt} = v_z \ , @f]
 * @f[ \frac{dv_z}{dt} = - \frac{e}{m_e} E_z(z,t) \left( 1- \frac{v_z^2}{c^2} \right)^{3/2} \ , @f]
 * with \f$ z(0) = 0 \f$ , \f$ v_z(0) = 0 \f$ .
 *
 * Since \f$ E_r = H_\phi = 0 \f$ at \f$ r=0 \f$ , the longitudinal electric field is directly computed inside this function from
   the expression
 * @f[ E_z =  \textrm{Re} \left\{ -2 E_0 e^{1/2} \left( \frac{a}{k_0} \right)^{3/2} \frac{1}{\tilde{R}} \left[
    \frac{1}{2} \left( 3\cos^2 \tilde{\theta} - 1 \right) \left(
    \frac{G_-^{(0)}}{\tilde{R}^2} - \frac{G_+^{(1)}}{c\tilde{R}} +  \frac{G_-^{(2)}}{3c^2} \right)  -
    \frac{G_-^{(2)}}{3c^2}\right] \right\}  \ .@f]
 * where \f$ \tilde{R} = z+ja \f$.
 *
 * @param t Current time \f$ t \f$
 * @param y[] Pointer to the position vector \f$  y=(z,v_z) \f$
 * @param f[] Pointer to a container that will hold \f$ dy/dt \f$
 * @param param Void pointer to the parameters. Need to be dereferenced into the appropriate structure
 * @return GSL_SUCCESS
 */
int dydt_npTM01_onaxis(double t, const double y[], double f[], void * param) {

  // Cast parameters
  NPparams& p = *static_cast<NPparams* >(param);

  // Create multi_array references
  typedef boost::multi_array_ref<const double,1> CSTDBLvecref;
  typedef boost::multi_array_ref<double,1> DBLvecref;
  typedef CSTDBLvecref::index index;
  CSTDBLvecref::extent_gen extents1;
  DBLvecref::extent_gen extents2;
  CSTDBLvecref Y(y,extents1[2]);
  DBLvecref F(f,extents2[2]);

  // Compute the longitudinal electric field component
  double psi0 = -exp(0.5)*pow(p.a*C/p.omega0,1.5)*p.E0;
  complex<double> ztilde = Y[0] + I*p.a;
  complex<double> tplus = t + ztilde/C + I*p.a/C;
  complex<double> tminus = t - ztilde/C + I*p.a/C;
  double fplusr = abs(1.0-I*p.omega0*tplus/p.s);
  double fplusphi = arg(1.0-I*p.omega0*tplus/p.s);
  double fminusr = abs(1.0-I*p.omega0*tminus/p.s);
  double fminusphi = arg(1.0-I*p.omega0*tminus/p.s);
  complex<double> Gminus0 = p.G0[0]*( pow(fplusr,-(p.s+1.0))*exp(-I*fplusphi*(p.s+1.0)) - pow(fminusr,-(p.s+1.0))*exp(-I*fminusphi*(p.s+1.0)) );
  complex<double> Gplus1 = p.G0[1]*( pow(fplusr,-(p.s+2.0))*exp(-I*fplusphi*(p.s+2.0)) + pow(fminusr,-(p.s+2.0))*exp(-I*fminusphi*(p.s+2.0)) );
  double Ez = real( (2.0*psi0/ztilde)*( Gminus0/pow(ztilde,2.0) - Gplus1/(C*ztilde) ) );

  // Equations of motion:
  F[0] = Y[1]; // dzdt
  F[1] = -(QE/ME)*Ez*pow(1-pow(Y[1],2.0)/pow(C,2.0),1.5); // dvz/dt

  return GSL_SUCCESS;

} // end function dydt_npTM01_onaxis

#endif // NPTM01_HPP_INCLUDED
