// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to calibrate the frequency of a simple
           sine wave based on noisy data.
*/

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Algorithm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <fstream>
#include <iostream>
#include <complex>
#include <valarray>

typedef double RealT;

typedef std::complex<RealT> Complex;
typedef std::valarray<Complex> CArray;
typedef std::valarray<RealT> RArray;


class CalibrationObjective : public ROL::Objective<RealT> {
private:
  const ROL::Ptr<std::vector<RealT> > data_; // vector of "measurements"
  const ROL::Ptr<std::vector<RealT> > time_; // time vector
  const RealT phase_;                            // wave phase
  const RealT amplitude_;                        // wave amplitude
  const RealT exp_const_;                        // exponential decay constant
  RealT data_scaling_;                           // scales the objective function

public:
  // Constructor.
  CalibrationObjective( const ROL::Ptr<std::vector<RealT> > & data,
                        const ROL::Ptr<std::vector<RealT> > & time,
                        RealT phase,
                        RealT amplitude,
                        RealT exp_const ) :
    data_(data), time_(time), phase_(phase), amplitude_(amplitude), exp_const_(exp_const) {

    unsigned num_samples = data_->size();
    data_scaling_ = 0.0;
    for (unsigned i=0; i<num_samples; ++i) {
      data_scaling_ += pow((*data_)[i], 2); 
    }
    if (data_scaling_ == 0.0) {
      data_scaling_ = 1.0;
    }

  }

  // Value of the calibration objective.
  RealT value( const ROL::Vector<RealT> &omega, RealT &tol ) {
    ROL::Ptr<const std::vector<RealT> > omega_vec_ptr =
      (dynamic_cast<const ROL::StdVector<RealT>&>(omega)).getVector();

    unsigned num_samples = data_->size();
    RealT y(0);
    RealT t(0);
    RealT val(0);
    RealT omega_scalar = (*omega_vec_ptr)[0];
    for (unsigned i=0; i<num_samples; ++i) {
      t   =  (*time_)[i];
      y   =  amplitude_*exp(-exp_const_*t)*sin(t*omega_scalar+phase_);
      val += pow(y-(*data_)[i], 2); 
    }
    return 0.5*val/data_scaling_;
  }

  // Gradient of the calibration objective.
  void gradient(ROL::Vector<RealT> &g, const ROL::Vector<RealT> &omega, RealT &tol ) {
    ROL::Ptr<std::vector<RealT> > g_vec_ptr =
      (dynamic_cast<ROL::StdVector<RealT>&>(g)).getVector();
    ROL::Ptr<const std::vector<RealT> > omega_vec_ptr =
      (dynamic_cast<const ROL::StdVector<RealT>&>(omega)).getVector();

    unsigned num_samples = data_->size();
    RealT gy(0);
    RealT t(0);
    RealT val(0);
    RealT omega_scalar = (*omega_vec_ptr)[0];
    for (unsigned i=0; i<num_samples; ++i) {
      t   =  (*time_)[i];
      gy  =  amplitude_*exp(-exp_const_*t)*t*cos(t*omega_scalar+phase_)*(amplitude_*exp(-exp_const_*t)*sin(t*omega_scalar+phase_) - (*data_)[i]);
      val += gy; 
    }
    (*g_vec_ptr)[0] = val/data_scaling_;
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
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

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

    // Our FFT is only robust for integer periods and numbers of samples that are powers of 2.
    int num_samples=pow(2,10);
    //RealT num_periods=12.3456;
    int num_periods=12;
    RealT epsilon = 3;
    RealT kappa=1;
    RealT phase = M_PI/7;
    RealT amplitude = 9.876;

    ROL::Ptr<std::vector<RealT> > data_ptr = ROL::makePtr<std::vector<RealT>>(num_samples, 0.0);
    ROL::Ptr<std::vector<RealT> > time_ptr = ROL::makePtr<std::vector<RealT>>(num_samples, 0.0);

    // This is for a decay
    RealT decay=200; // number of periods for an e decay
    RealT noise_frequency = 2*4.5; // ratio of this to real frequency
    RealT noise_amplitude = 12*0.124; // ratio of this to real amplitude
    //RealT noise_frequency(0); // ratio of this to real frequency
    //RealT noise_amplitude(0); // ratio of this to real amplitude
    RealT noise_phase = M_PI*exp(1.);

    RealT epsilon_0=8.854187817e-12;
    RealT mu=1.2566370614e-6;
    RealT c = 1./sqrt(epsilon*epsilon_0*mu);

    // omega*kappa = c
    RealT omega = c/kappa;
    RealT frequency = omega/(2*M_PI);
    RealT period = 1./frequency;
    RealT total_time = num_periods*period;
    RealT dt = total_time/(num_samples-1);
    RealT exp_const = 1/(decay*period);
    RealT noise_omega = noise_frequency*omega;

    // Generate "measurements" and output to file.
    *outStream << std::endl << "Generating measurements through a noisy computation:" << std::endl;
    CArray data_fft(num_samples);
    std::ofstream measfile;
    measfile.open("measurements.txt");
    for (int i=0; i<num_samples; ++i) {
      RealT t = dt*i;
      //RealT t = total_time*(RealT)rand()/(RealT)RAND_MAX;
      (*time_ptr)[i] = t;
      // additive noise
      (*data_ptr)[i] = amplitude*(sin(t*omega+phase)+ noise_amplitude*sin(t*noise_omega+noise_phase))*exp(-exp_const*t);
      (*data_ptr)[i] += amplitude*(0.5*noise_amplitude*sin(t*0.5*noise_omega+noise_phase))*exp(-exp_const*t);
      (*data_ptr)[i] += amplitude*(0.4*noise_amplitude*sin(t*0.4*noise_omega+noise_phase))*exp(-exp_const*t);
      (*data_ptr)[i] += amplitude*(0.3*noise_amplitude*sin(t*0.3*noise_omega+noise_phase))*exp(-exp_const*t);
      (*data_ptr)[i] += amplitude*(0.2*noise_amplitude*sin(t*0.2*noise_omega+noise_phase))*exp(-exp_const*t);
      (*data_ptr)[i] += amplitude*(0.1*noise_amplitude*sin(t*0.1*noise_omega+noise_phase))*exp(-exp_const*t);
      measfile   << "   " << std::scientific << std::left << std::setprecision(4) << std::setw(6) << t 
                 << "   " << std::right << std::setw(14) << (*data_ptr)[i] << std::endl;
      data_fft[i].real((*data_ptr)[i]);
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
      omega_i = 2*M_PI*i/((num_samples-1)*dt);
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
    ROL::Ptr<std::vector<RealT> > omega_vec_ptr = ROL::makePtr<std::vector<RealT>>(1, 0.0);
    ROL::StdVector<RealT> omega_rol_vec(omega_vec_ptr);

    // Define calibration objective.
    CalibrationObjective cal_obj(data_ptr, time_ptr, phase, amplitude, exp_const);

    // Output objective values to file given a sweep of parameters.
    *outStream << std::endl << "Sampling calibration objective function:" << std::endl;
    RealT omega_first = 0.0*omega;
    RealT omega_last = 10*omega;
    int num_omegas = 1000;
    RealT tol = 0.0;
    std::ofstream objfile;
    objfile.open("cal_obj.txt");
    for (int i=0; i<num_omegas; ++i) {
      RealT omega_step = (omega_last-omega_first)/(num_omegas-1)*i+omega_first;
      (*omega_vec_ptr)[0] = omega_step;
      objfile << "   " <<  std::scientific << std::left << std::setprecision(4) << std::setw(6) << omega_step 
              << "   " << std::right << std::setw(14) << cal_obj.value(omega_rol_vec, tol) << std::endl;      
    }
    *outStream << "   Done.  You can visualize results in gnuplot:    plot \"cal_obj.txt\" with lines" << std::endl;
    objfile.close();

    *outStream << std::endl << "Matching measurements to calibrate frequency:" << std::endl;
    
    // Define ROL step, status test, and algorithm.
    RealT gtol     = 1e-16;          // gradient tolerance
    RealT stol     = 1e-18;          // step tolerance
    int   max_iter = 100;            // maximum number of optimization iterations
    ROL::ParameterList parlist;  // list of algorithmic parameters
      parlist.sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type", "Quasi-Newton Method");
      parlist.sublist("General").sublist("Secant").set("Type", "Limited-Memory BFGS");
      // Trust-region step parameters.
      parlist.sublist("Step").sublist("Trust Region").set("Subproblem Solver", "Truncated CG");
      parlist.sublist("Step").sublist("Trust Region").set("Initial Radius", 1e7);
      parlist.sublist("Step").sublist("Trust Region").set("Maximum Radius", 1e12);
    ROL::Ptr<ROL::LineSearchStep<RealT> >   lsstep = ROL::makePtr<ROL::LineSearchStep<RealT>>(parlist);  // line-search method
    ROL::Ptr<ROL::TrustRegionStep<RealT> >  trstep = ROL::makePtr<ROL::TrustRegionStep<RealT>>(parlist);  // trust-region method
    ROL::Ptr<ROL::StatusTest<RealT> >       status = ROL::makePtr<ROL::StatusTest<RealT>>(gtol, stol, max_iter);  // status test

    // Run simple algorithm (starting at many initial points).
    /*
    int   num_inits     = 1000;
    RealT omega_min     = 1e6;
    RealT omega_max     = 1e10;
    RealT objval_min    = std::numeric_limits<RealT>::max();
    RealT solution_min  = std::numeric_limits<RealT>::max();
    for (int i=0; i<num_inits; ++i) {  // start at several initial guesses
      (*omega_vec_ptr)[0] = omega_min + (RealT)rand() / ((RealT)RAND_MAX/(omega_max-omega_min)); 
      ROL::Algorithm<RealT> algo(trstep, status, false);
      algo.run(omega_rol_vec, cal_obj, true, *outStream);
      RealT solution = (*omega_vec_ptr)[0];
      RealT objval   = cal_obj.value(omega_rol_vec, tol);
      if (objval < objval_min) {
        objval_min = objval;
        solution_min = solution;
      } 
    }
    */


    /*** Multilevel single linkage, where all local minima are found
         almost surely, at the expense of a finite number of local
         optimization runs, almost surely. ***/ 

    RealT omega_min = 0.0*omega;
    RealT omega_max = 10.0*omega;
    RealT realmax = std::numeric_limits<RealT>::max();
    RealT objval_min = realmax;
    RealT solution_min = realmax;
    int dim = 1;
    int k_max = 20;
    int num_points = 150;
    RealT r_k(0);
    RealT sigma(4.0);
    RealT dist_to_loc(10*(omega_max-omega_min)/(k_max*num_points));
    srand(0);
    std::vector<ROL::Ptr<ROL::Vector<RealT > > >  vec_sample;
    std::vector<ROL::Ptr<ROL::Vector<RealT > > >  vec_locmin;
    std::vector<RealT> val_sample;
    std::vector<RealT> val_locmin;
    std::vector<RealT> min_distance;
    ROL::Ptr<ROL::Vector<RealT> > tmp_vec  = omega_rol_vec.clone();

    for (int k=0; k<k_max; ++k) {

      for (int i=0; i<num_points; ++i) {
        vec_sample.push_back(omega_rol_vec.clone());
        // Compute random sample ... this would have to be generalized.
        //(vec_sample.back())->randomize();
          RealT tmp = omega_min + (RealT)rand() / ((RealT)RAND_MAX/(omega_max-omega_min));
          ROL::Ptr<std::vector<RealT> > last_vec_ptr =
            (dynamic_cast<ROL::StdVector<RealT>&>(*(vec_sample.back()))).getVector();
          (*last_vec_ptr)[0] = tmp;

        // Compute objective function value at the sample.
        val_sample.push_back(cal_obj.value(*(vec_sample.back()), tol));

        // Initialize minimum distance to other points to numeric max.
        min_distance.push_back(realmax);

        // Compute minimum distances to points in the sample set.
        std::vector<RealT> tmp_distance;
        for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itsam = vec_sample.begin(); itsam != vec_sample.end(); ++itsam) {
          tmp_vec->set(*(vec_sample.back()));
          tmp_vec->axpy(-1.0, **itsam);
          RealT dist  = tmp_vec->norm();
          tmp_distance.push_back(dist);
        }
        tmp_distance.back() = realmax;
        for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itsam = vec_sample.begin(); itsam != vec_sample.end(); ++itsam) {
          int idx = itsam - vec_sample.begin(); // current iterator index
          if ((tmp_distance[idx] < min_distance[idx]) && (itsam != vec_sample.end()-1)) {
            min_distance[idx] = tmp_distance[idx];
          }
        }
        min_distance.back() = *std::min_element(tmp_distance.begin(),tmp_distance.end());
      }

      // Compute current minimum over all samples.
      RealT minval_sample = *std::min_element(val_sample.begin(),val_sample.end());

      // Compute r_k based on the domain volume ... this would have to be generalized.
      int nsamples = vec_sample.size();
      //r_k = (1.0/M_PI)*pow(std::tgamma(1.0+dim/2.0)*volumeD*sigma*log10(nsamples)/nsamples, 1.0/dim);
      //r_k = (1.0/sqrt(M_PI))*pow(std::tgamma(1.0+dim/2.0)*(omega_max-omega_min)*sigma*log10(nsamples)/nsamples, 1.0/dim);
        r_k = (1.0/sqrt(M_PI))*pow(tgamma(1.0+dim/2.0)*(omega_max-omega_min)*sigma*log10(nsamples)/nsamples, 1.0/dim);

      // Start local optimization runs.
      for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itsam = vec_sample.begin(); itsam != vec_sample.end(); ++itsam) {
        bool islocal = false;
        bool isnearlocal = false;
        for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itloc = vec_locmin.begin(); itloc != vec_locmin.end(); ++itloc) {
          if (*itsam == *itloc) {
            islocal = true;
          }
        }
        for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itloc = vec_locmin.begin(); itloc != vec_locmin.end(); ++itloc) {
          tmp_vec->set(**itloc);
          tmp_vec->axpy(-1.0, **itsam);
          RealT dist  = tmp_vec->norm();
          if (dist < dist_to_loc) {
              isnearlocal = true;
          }
        }
        if (!islocal && !isnearlocal) {
          int idx = itsam - vec_sample.begin(); // current iterator index
          if ((val_sample[idx] <= minval_sample) || (min_distance[idx] > r_k)) {
            ROL::Algorithm<RealT> algo(trstep, status, false);
            ROL::Ptr<ROL::Vector<RealT> > soln_vec  = omega_rol_vec.clone();
            soln_vec->set(**itsam);
            algo.run(*soln_vec, cal_obj, true, *outStream);
            vec_locmin.push_back(*itsam);
            val_locmin.push_back(cal_obj.value(*soln_vec, tol));
          }
        }
      } // end for-loop over samples
    }  // end outer iteration k

    *outStream << std::endl << "Number of local minima identified: " << val_locmin.size() << std::endl;
    *outStream << "Minimizers:" << std::endl;
    for (std::vector<ROL::Ptr<ROL::Vector<RealT > > >::iterator itloc = vec_locmin.begin(); itloc != vec_locmin.end(); ++itloc) {
      ROL::Ptr<const std::vector<RealT> > vec_ptr =
        (dynamic_cast<ROL::StdVector<RealT>&>(**itloc)).getVector();
      *outStream << "  " << (*vec_ptr)[0] << std::endl;
    }
    *outStream << std::endl;

    // Extract minimum and minimizer.  
    std::vector<RealT>::iterator itmin = std::min_element(val_locmin.begin(),val_locmin.end());
    objval_min = *itmin;
    int idx_min = itmin - val_locmin.begin();
    solution_min = (*(dynamic_cast<ROL::StdVector<RealT>&>(*(vec_locmin[idx_min]))).getVector())[0];

    /*** End of MLSL. ***/

    *outStream << std::endl << "'True' frequency:             omega = " << std::left << std::scientific << std::setprecision(8) << std::setw(12) << omega; 
    *outStream << std::endl << "FFT frequency:            omega_fft = " << maxomega; 
    *outStream << std::endl << "Computed optimal frequency:  omega* = " << solution_min << std::endl; 
    *outStream << std::endl << "'True' epsilon:                 eps = " << std::left << std::setprecision(8) << std::setw(12) << epsilon; 
    *outStream << std::endl << "FFT epsilon:                eps_fft = " << 1/(pow(maxomega*kappa,2)*epsilon_0*mu); 
    *outStream << std::endl << "Computed optimal epsilon       eps* = " << 1/(pow(solution_min*kappa,2)*epsilon_0*mu) << std::endl; 
    (*omega_vec_ptr)[0] = maxomega;
    *outStream << std::endl << "Objective value at FFT freq     val = " << cal_obj.value(omega_rol_vec, tol); 
    *outStream << std::endl << "Objective value at opt freq    val* = " << objval_min << std::endl; 

    if ( std::abs(solution_min - omega)/omega > 1e-2 ) {
      *outStream << std::endl << "WARNING: Relative error = " << std::abs(solution_min - omega)/omega << std::endl; 
      errorFlag += 1;
    }

  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;

}

