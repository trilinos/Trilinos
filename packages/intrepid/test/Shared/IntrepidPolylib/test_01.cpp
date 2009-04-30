// @HEADER
///////////////////////////////////////////////////////////////////////////////
//
// File: test_01.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
// This file is redistributed with the Intrepid package. It should be used
// in accordance with the above MIT license, at the request of the authors.
// This file is NOT covered by the usual Intrepid/Trilinos LGPL license.
//
// Origin: Nektar++ library, http://www.nektar.info, downloaded on
//         March 10, 2009.
//
///////////////////////////////////////////////////////////////////////////////


/** \file   test_01.cpp
    \brief  Test file for a set of functions providing orthogonal polynomial
            polynomial calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
            modified and redistributed by D. Ridzal.
*/

#include "Intrepid_Polylib.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"


using namespace std;
using namespace Intrepid;

#define NPLOWER  5
#define NPUPPER 15
#define TEST_EPS  1000*INTREPID_TOL

#define GAUSS_INT            1
#define GAUSS_RADAUM_INT     1
#define GAUSS_RADAUP_INT     1
#define GAUSS_LOBATTO_INT    1
#define GAUSS_DIFF           1
#define GAUSS_RADAUM_DIFF    1
#define GAUSS_RADAUP_DIFF    1
#define GAUSS_LOBATTO_DIFF   1
#define GAUSS_INTERP         1
#define GAUSS_RADAUM_INTERP  1
#define GAUSS_RADAUP_INTERP  1
#define GAUSS_LOBATTO_INTERP 1


#define INTREPID_TEST_COMMAND( S )                                                                                  \
{                                                                                                                   \
  try {                                                                                                             \
    S ;                                                                                                             \
  }                                                                                                                 \
  catch (std::logic_error err) {                                                                                    \
      *outStream << "Expected Error ----------------------------------------------------------------\n";            \
      *outStream << err.what() << '\n';                                                                             \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";    \
  };                                                                                                                \
}

/* local routines */
template<class Scalar>
Scalar ddot(int n, Scalar *x, int incx, Scalar *y, int incy)
{
  register Scalar sum = 0.;

  while (n--) {
    sum += (*x) * (*y);
    x   += incx;
    y   += incy;
  }
  return sum;
}

template<class Scalar>
Scalar * dvector(int nl,int nh)
{
  Scalar * v;

  v = (Scalar *)malloc((unsigned) (nh-nl+1)*sizeof(Scalar));
  return v-nl;
}


/* -------------------------------------------------------------------
   This is a routine to test the integration, differentiation and
   interpolation routines in the polylib.c.

   First, it performs the integral

      /1      alpha   beta  alpha,beta
     |   (1-x)   (1+x)     P (x)       dx  = 0
     /-1                    n

   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= 2*np - delta

   delta = 1 (gauss), 2(radau), 3(lobatto).
   The integral is evaluated and if it is larger that EPS then the
   value of alpha,beta,np,n and the integral is printed to the screen.

   After every alpha value the statement
       "finished checking all beta values for alpha = #"
   is printed

   The routine then evaluates the derivate of

          d   n      n-1
      -- x  = n x
      dx

   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= np - 1

   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
       "finished checking all beta values for alpha = #"
   is printed

   Finally the routine  evaluates the interpolation of

             n      n
        z  to  x

   where z are the quadrature zeros and x are the equispaced points

                  2*i
        x    =   -----   - 1.0    (0 <= i <= np-1)
     i       (np-1)


   for all   -0.5 <= alpha <= 5   (increments of 0.5)
             -0.5 <= beta  <= 5   (increments of 0.5)

   using np points where
          NPLOWER <= np <= NPUPPER
                2     <= n  <= np - 1

   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
      "finished checking all beta values for alpha = #"
   is printed

   The above checks are performed for all the Gauss, Gauss-Radau and
   Gauss-Lobatto points. If you want to disable any routine then set
      GAUSS_INT, GAUSS_RADAU_INT, GAUSS_LOBATTO_INT = 0
   for the integration rouintes
      GAUSS_DIFF,GAUSS_RADAU_DIFF, GAUSS_LOBATTO_DIFF = 0
   for the differentiation routines
      GAUSS_INTERP,GAUSS_RADAU_INTERP, GAUSS_LOBATTO_INTERP = 0
   for the interpolation routines.
------------------------------------------------------------------*/
int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Unit Test (IntrepidPolylib)                           |\n" \
  << "|                                                                             |\n" \
  << "|     1) the original Polylib tests                                           |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  int errorFlag = 0;
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 1;

  typedef IntrepidPolylib ipl; 
  IntrepidPolylib iplib;

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: exceptions                                                          |\n"\
  << "===============================================================================\n";

  try{
    INTREPID_TEST_COMMAND( iplib.gammaF((double)0.33) );
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };

  if (TestForException_getThrowNumber() != endThrowNumber)
    errorFlag = -1000;

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try {
      { // start scope
      int np,n,i;
      double *z,*w,*p,sum=0,alpha,beta,*d;

      z  = dvector<double>(0,NPUPPER-1);
      w  = dvector<double>(0,NPUPPER-1);
      p  = dvector<double>(0,NPUPPER-1);

      d  = dvector<double>(0,NPUPPER*NPUPPER-1);

#if GAUSS_INT
      /* Gauss Integration */
      *outStream << "Begin checking Gauss integration\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgj(z,w,np,alpha,beta);
            for(n = 2; n < 2*np-1; ++n){
              ipl::jacobfd(np,z,p,(double*)0,n,alpha,beta);
              sum = ddot(np,w,1,p,1);
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  integral was " << sum << "\n";
              }
            }
          }

          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss Integration\n";
#endif

#if GAUSS_RADAUM_INT
      /* Gauss Radau (z = -1) Integration */
      *outStream << "Begin checking Gauss-Radau (z = -1) integration\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjm(z,w,np,alpha,beta);
            for(n = 2; n < 2*np-2; ++n){
              ipl::jacobfd(np,z,p,(double*)0,n,alpha,beta);
              sum = ddot(np,w,1,p,1);
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  integral was " << sum << "\n";
              }
            }
          }

          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z = -1) Integration\n";
#endif

#if GAUSS_RADAUP_INT
      /* Gauss Radau (z = +1) Integration */
      *outStream << "Begin checking Gauss-Radau (z = +1) integration\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjp(z,w,np,alpha,beta);
            for(n = 2; n < 2*np-2; ++n){
              ipl::jacobfd(np,z,p,(double*)0,n,alpha,beta);
              sum = ddot(np,w,1,p,1);
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  integral was " << sum << "\n";
              }
            }
          }

          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z = +1) Integration\n";
#endif

#if GAUSS_LOBATTO_INT
      /* Gauss Lobatto Integration */
      *outStream << "Begin checking Gauss-Lobatto integration\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwglj(z,w,np,alpha,beta);
            for(n = 2; n < 2*np-3; ++n){
              ipl::jacobfd(np,z,p,(double*)0,n,alpha,beta);
              sum = ddot(np,w,1,p,1);
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  integral was " << sum << "\n";
              }
            }
          }

          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Lobatto Integration\n";
#endif

#if GAUSS_DIFF
      *outStream << "Begin checking differentiation through Gauss points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgj(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              ipl::Dgj(d,z,np,alpha,beta);

              for(i = 0; i < np; ++i) p[i] = std::pow(z[i],n);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - n*std::pow(z[i],n-1));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss Jacobi differentiation\n";
#endif

#if GAUSS_RADAUM_DIFF
      *outStream << "Begin checking differentiation through Gauss-Radau (z=-1) points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjm(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              ipl::Dgrjm(d,z,np,alpha,beta);

              for(i = 0; i < np; ++i) p[i] = std::pow(z[i],n);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - n*std::pow(z[i],n-1));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z=-1) differentiation\n";
#endif

#if GAUSS_RADAUP_DIFF
      *outStream << "Begin checking differentiation through Gauss-Radau (z=+1) points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjp(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              ipl::Dgrjp(d,z,np,alpha,beta);

              for(i = 0; i < np; ++i) p[i] = std::pow(z[i],n);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - n*std::pow(z[i],n-1));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z=+1) differentiation\n";
#endif

#if GAUSS_RADAUP_DIFF
      *outStream << "Begin checking differentiation through Gauss-Lobatto (z=+1) points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwglj(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              ipl::Dglj(d,z,np,alpha,beta);

              for(i = 0; i < np; ++i) p[i] = std::pow(z[i],n);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - n*std::pow(z[i],n-1));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Lobatto differentiation\n";
#endif

#if GAUSS_INTERP
      *outStream << "Begin checking interpolation through Gauss points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgj(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              for(i = 0; i < np; ++i) {
                w[i] = 2.0*i/(double)(np-1)-1.0;
                p[i] = std::pow(z[i],n);
              }
              ipl::Imgj(d,z,w,np,np,alpha,beta);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - std::pow(w[i],n));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss Jacobi interpolation\n";
#endif

#if GAUSS_RADAUM_INTERP
      *outStream << "Begin checking interpolation through Gauss-Radau (z=-1) points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjm(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              for(i = 0; i < np; ++i) {
                w[i] = 2.0*i/(double)(np-1)-1.0;
                p[i] = std::pow(z[i],n);
              }
              ipl::Imgrjm(d,z,w,np,np,alpha,beta);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - std::pow(w[i],n));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z=-1) interpolation\n";
#endif

#if GAUSS_RADAUP_INTERP
      *outStream << "Begin checking interpolation through Gauss-Radau (z=+1) points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwgrjp(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              for(i = 0; i < np; ++i) {
                w[i] = 2.0*i/(double)(np-1)-1.0;
                p[i] = std::pow(z[i],n);
              }
              ipl::Imgrjp(d,z,w,np,np,alpha,beta);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - std::pow(w[i],n));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Radau (z=+1) interpolation\n";
#endif

#if GAUSS_LOBATTO_INTERP
      *outStream << "Begin checking interpolation through Gauss-Lobatto points\n";
      alpha = -0.5;
      while(alpha <= 5.0){
        beta = -0.5;
        while(beta <= 5.0){

          for(np = NPLOWER; np <= NPUPPER; ++np){
            ipl::zwglj(z,w,np,alpha,beta);
            for(n = 2; n < np-1; ++n){
              for(i = 0; i < np; ++i) {
                w[i] = 2.0*i/(double)(np-1)-1.0;
                p[i] = std::pow(z[i],n);
              }
              ipl::Imglj(d,z,w,np,np,alpha,beta);
              sum = 0;
              for(i = 0; i < np; ++i)
                sum += std::abs(ddot(np,d+i*np,1,p,1) - std::pow(w[i],n));
              sum /= np;
              if(std::abs(sum)>TEST_EPS) {
                errorFlag = -1000;
                *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                              ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
              }
            }
          }
          beta += 0.5;
        }
        *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
        alpha += 0.5;
      }
      *outStream << "Finished checking Gauss-Lobatto interpolation\n";
#endif

      free(z); free(w); free(p); free(d);

      } // end scope

      /******************************************/
      *outStream << "\n";
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };


  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
