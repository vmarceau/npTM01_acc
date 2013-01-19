/**
 * @file    npTM01_readdata.cpp
 * @brief   Source file used to read data files

 * This source file is used to read and analyze data files when then are too large for these tasks to be performed in standard
 * softwares such as Matlab and GNU Octave. The code is written in order to minimize the amount of memory required, and can
 * therefore require a significant computation time.
 * @author  Vincent Marceau (vincent.marceau.2@ulaval.ca)
 * @since   January 2012
 * @date    January 2012
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

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
int main() {

  int numheaderlines = 2;
  int numrows = 3;

  ifstream datafile;
  datafile.open ("./dat/readdata/datatest.dat", ios::in);

  for ( int n=0; n<numheaderlines; n++ ) {
    string dummyline;
    getline (datafile,dummyline);
  } // end for

  double z,phi,E;
  double Emax = 0;
  double Emin = 0;
  while ( datafile.good() ) {

    for (int m=0; m<numrows; m++) {
      datafile >> z;
      datafile >> phi;
      datafile >> E;
      if ( !datafile.good() )
        break;
      if ( E > Emax )
        Emax = E;
      if ( E < Emin )
        Emin = E;

      cout << "Emax = " << Emax << endl;
      cout << "Emin = " << Emin << endl;

    } // end for

  } // end while

  cout << "Emax = " << Emax << endl;
  cout << "Emin = " << Emin << endl;

  return 0;

} // end function main
