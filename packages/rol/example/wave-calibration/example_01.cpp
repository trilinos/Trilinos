// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to calibrate the frequency of a simple
           sine wave based on noisy data.
*/

#include "ROL_Objective.hpp"
#include "ROL_EqualityConstraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <complex>
#include <valarray>

typedef double RealT;

typedef std::complex<RealT> Complex;
typedef std::valarray<Complex> CArray;
typedef std::valarray<RealT> RArray;


class CalibrationObjective : public ROL::Objective<RealT> {
private:
  const Teuchos::RCP<std::vector<RealT> > data_; // vector of "measurements"
  const Teuchos::RCP<std::vector<RealT> > time_; // time vector
  const RealT phase_;                            // wave phase
  const RealT amplitude_;                        // wave amplitude
  const RealT exp_const_;                        // exponential decay constant
  RealT data_scaling_;                           // scales the objective function

public:
  // Constructor.
  CalibrationObjective( const Teuchos::RCP<std::vector<RealT> > & data,
                        const Teuchos::RCP<std::vector<RealT> > & time,
                        RealT phase,
                        RealT amplitude,
                        RealT exp_const ) :
    data_(data), time_(time), phase_(phase), amplitude_(amplitude), exp_const_(exp_const) {

    unsigned num_samples = data_->size();
    data_scaling_ = 0.0;
    for (unsigned i=0; i<num_samples; ++i) {
      data_scaling_ += pow((*data_)[i], 2); 
    }

  }

  // Value of the calibration objective.
  RealT value( const ROL::Vector<RealT> &omega, RealT &tol ) {
    Teuchos::RCP<const std::vector<RealT> > omega_vec_rcp =
      (Teuchos::dyn_cast<const ROL::StdVector<RealT> >(omega)).getVector();

    unsigned num_samples = data_->size();
    RealT y(0);
    RealT t(0);
    RealT val(0);
    RealT omega_scalar = (*omega_vec_rcp)[0];
    for (unsigned i=0; i<num_samples; ++i) {
      t   =  (*time_)[i];
      y   =  amplitude_*exp(-exp_const_*t)*sin(t*omega_scalar+phase_);
      val += pow(y-(*data_)[i], 2); 
    }
    return 0.5*val/data_scaling_;
  }

  // Gradient of the calibration objective.
  void gradient(ROL::Vector<RealT> &g, const ROL::Vector<RealT> &omega, RealT &tol ) {
    Teuchos::RCP<std::vector<RealT> > g_vec_rcp =
      (Teuchos::dyn_cast<ROL::StdVector<RealT> >(g)).getVector();
    Teuchos::RCP<const std::vector<RealT> > omega_vec_rcp =
      (Teuchos::dyn_cast<const ROL::StdVector<RealT> >(omega)).getVector();

    unsigned num_samples = data_->size();
    RealT gy(0);
    RealT t(0);
    RealT val(0);
    RealT omega_scalar = (*omega_vec_rcp)[0];
    for (unsigned i=0; i<num_samples; ++i) {
      t   =  (*time_)[i];
      gy  =  amplitude_*exp(-exp_const_*t)*t*cos(t*omega_scalar+phase_)*(amplitude_*exp(-exp_const_*t)*sin(t*omega_scalar+phase_) - (*data_)[i]);
      val += gy; 
    }
    (*g_vec_rcp)[0] = val/data_scaling_;
  }

};

// Cooley–Tukey FFT (in-place, divide-and-conquer).
// Higher memory requirements and redundancy although more intuitive.
// Taken from Wikipedia.
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;
 
    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];
 
    // conquer
    fft(even);
    fft(odd);
 
    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  // *** Example body.

  try {

    *outStream << std::endl << std::endl << 
    "Given the equation" << std::endl << std::endl <<
    "  y(t) = A*sin(omega*t+phi)*exp(-D*t)" << std::endl << std::endl <<
    "where A is the amplitude, omega is the angular frequency," << std::endl <<
    "phi is the phase and D is the decay constant, we calibrate" << std::endl <<
    "the angular frequency omega based on noisy measurements" << std::endl <<
    "y_meas(t) of y(t) at discrete points in time." << std::endl;

    int num_samples=pow(2,8);
    double num_periods=12.3456;
    double epsilon = 3;
    double k=1;
    double phase = M_PI/7;
    double amplitude = 9.876;

    Teuchos::RCP<std::vector<RealT> > data_rcp = Teuchos::rcp( new std::vector<RealT> (num_samples, 0.0) );
    Teuchos::RCP<std::vector<RealT> > time_rcp = Teuchos::rcp( new std::vector<RealT> (num_samples, 0.0) );

    // This is for a decay
    double decay=200; // number of periods for an e decay
    double noise_frequency = 4.5; // ratio of this to real frequency
    double noise_amplitude = 0.124; // ratio of this to real amplitude
    //double noise_frequency = 0; // ratio of this to real frequency
    //double noise_amplitude = 0; // ratio of this to real amplitude
    double noise_phase = M_PI*exp(1.);

    double epsilon_0=8.854187817e-12;
    double mu=1.2566370614e-6;
    double c = 1./sqrt(epsilon*epsilon_0*mu);

    // omega k = c
    double omega = c/k;
    double frequency = omega/(2*M_PI);
    double period = 1./frequency;
    double total_time = num_periods*period;
    double dt = total_time/(num_samples-1);
    double exp_const = 1/(decay*period);
    double noise_omega = noise_frequency*omega;

    // Generate "measurements" and output to file.
    *outStream << std::endl << "Generating measurements through a noisy computation:" << std::endl;
    CArray data_fft(num_samples);
    std::ofstream measfile;
    measfile.open("measurements.txt");
    for (int i=0; i<num_samples; ++i) {
      double t = dt*i;
      //double t = total_time*(double)rand()/(double)RAND_MAX;
      (*time_rcp)[i] = t;
      (*data_rcp)[i] = amplitude*(sin(t*omega+phase)+ noise_amplitude*sin(t*noise_omega+noise_phase))*exp(-exp_const*t);
      measfile   << "   " <<  std::scientific << std::left << std::setprecision(4) << std::setw(6) << t 
                 << "   " << std::right << std::setw(14) << (*data_rcp)[i] << std::endl;
      data_fft[i].real((*data_rcp)[i]);
      data_fft[i].imag(0);
    }
    *outStream << "   Done.  You can visualize results in gnuplot:    plot \"measurements.txt\" with lines" << std::endl;
    measfile.close();

    // Peform discrete Fourier transform, using an FFT algorithm, and output to file.
    // Also compute the frequency at which the maximum magnitude FFT is reached.
    *outStream << std::endl << "Performing FFT and computing frequency with max FFT magnitude:" << std::endl;
    fft(data_fft);
    data_fft /= num_samples;
    std::ofstream fftfile;
    fftfile.open("fft.txt");
    CArray magn(std::abs(data_fft));
    RArray magn_r(num_samples/2);
    RealT maxmagn(0);
    RealT maxomega(0);
    RealT omega_i(0);
    for (int i = 0; i < num_samples/2; ++i)
    {
      omega_i = 2*M_PI*i/(num_samples*dt);
      magn_r[i] = magn[i].real();
      // compute frequency with max FFT magnitude
      if (maxmagn < magn_r[i]) {
        maxmagn = magn_r[i];
        maxomega = omega_i;
      }
      fftfile   << "   " <<  std::scientific << std::left << std::setprecision(4) << std::setw(6) << omega_i
                << "   " << std::right << std::setw(14) << magn_r[i] << std::endl;
    }
    fftfile.close();
    *outStream << "   Done.  You can visualize results in gnuplot:    plot \"fft.txt\" with lines" << std::endl;
    fftfile.close();

    /************** Solve calibration problem. *************/

    // Define optimization 'vector' by using ROL::StdVector.
    Teuchos::RCP<std::vector<RealT> > omega_vec_rcp = Teuchos::rcp( new std::vector<RealT> (1, 0.0) );
    ROL::StdVector<RealT> omega_rol_vec(omega_vec_rcp);

    // Define calibration objective.
    CalibrationObjective cal_obj(data_rcp, time_rcp, phase, amplitude, exp_const);

    // Output objective values to file given a sweep of parameters.
    *outStream << std::endl << "Sampling calibration objective function:" << std::endl;
    double omega_first = 0.0*omega;
    double omega_last = 2*omega;
    int num_omegas = 10000;
    double tol = 0.0;
    std::ofstream objfile;
    objfile.open("cal_obj.txt");
    for (int i=0; i<num_omegas; ++i) {
      double omega_step = (omega_last-omega_first)/(num_omegas-1)*i+omega_first;
      (*omega_vec_rcp)[0] = omega_step;
      objfile << "   " <<  std::scientific << std::left << std::setprecision(4) << std::setw(6) << omega_step 
              << "   " << std::right << std::setw(14) << cal_obj.value(omega_rol_vec, tol) << std::endl;      
    }
    *outStream << "   Done.  You can visualize results in gnuplot:    plot \"cal_obj.txt\" with lines" << std::endl;
    objfile.close();

    *outStream << std::endl << "Matching measurements to calibrate frequency:" << std::endl;
    
    // Set initial guess for omega.
    RealT omega_init = maxomega;
    (*omega_vec_rcp)[0] = omega_init;

    // Define ROL step, status test, and algorithm.
    RealT gtol     = 1e-12;          // gradient tolerance
    RealT stol     = 1e-14;          // step tolerance
    int   max_iter = 100;            // maximum number of optimization iterations
    Teuchos::ParameterList parlist;  // list of algorithmic parameters
      // Line-search step parameters.
      parlist.set("Descent Type", "Quasi-Newton");
      parlist.set("Secant Type", "Limited-Memory BFGS");
      // Trust-region step parameters.
      parlist.set("Trust-Region Subproblem Solver Type", "Truncated CG");
      parlist.set("Initial Trust-Region Radius", 1e7);
      parlist.set("Maximum Trust-Region Radius", 1e12);
    ROL::LineSearchStep<RealT>   lsstep(parlist);  // line-search method
    ROL::TrustRegionStep<RealT>  trstep(parlist);  // trust-region method
    ROL::StatusTest<RealT>       status(gtol, stol, max_iter);
    ROL::DefaultAlgorithm<RealT> algo(trstep, status, false);

    // Run algorithm.
    algo.run(omega_rol_vec, cal_obj, true, *outStream);

    *outStream << std::endl << "'True' frequency:             omega = " << std::left << std::scientific << std::setprecision(8) << std::setw(12) << omega; 
    *outStream << std::endl << "FFT frequency:            omega_fft = " << maxomega; 
    *outStream << std::endl << "Initial frequency:       omega_init = " << omega_init; 
    *outStream << std::endl << "Computed optimal frequency:  omega* = " << (*omega_vec_rcp)[0] << std::endl; 
    *outStream << std::endl << "'True' epsilon:                 eps = " << std::left << std::setprecision(8) << std::setw(12) << epsilon; 
    *outStream << std::endl << "FFT epsilon:                eps_fft = " << 1/(pow(maxomega*k,2)*epsilon_0*mu); 
    *outStream << std::endl << "Computed optimal epsilon       eps* = " << 1/(pow((*omega_vec_rcp)[0]*k,2)*epsilon_0*mu) << std::endl; 

    if ( abs((*omega_vec_rcp)[0] - omega)/omega > 1e-3 ) {
      *outStream << std::endl << "WARNING: Relative error = " << abs((*omega_vec_rcp)[0] - omega)/omega << std::endl; 
      errorFlag += 1;
    }

  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

