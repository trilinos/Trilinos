// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "element_density_function.h"
#include <sstream>
#include <stdio.h>
namespace PAMGEN_NEVADA {


  /*****************************************************************************/
  void Element_Density_Function::Display_Class(
      std::ostream& s,
      const std::string &/* indent */
      )
  {
    s << std::endl
      << "Evaluation information for user defined run-time compiled " << std::endl
      << "Element Density Function" << std::endl
      << "For Direction " << application_direction << std::endl;
    s  << "Maximum function value = "
      << max_eval << std::endl << "Minimum function value = "
      << min_eval << std::endl
      << "Ratio  = " << max_eval/min_eval << std::endl
      << "Function integrates to " << integral_total/running_sum_length << " over mesh range." << std::endl << std::endl;
    double input_var = 0.;
    double * ivp = & input_var;
    _function.varAddrFill(0, ivp);
    double output_var = 0;
    double * ovp = & output_var;
    _function.varAddrFill(1, ovp);


#define Y_POINTS 20
#define X_POINTS 62

    //                   ict jct
    //                   x   y
    char gchars [X_POINTS][Y_POINTS];
    double delta_var = max_eval_range-min_eval_range;
    double gdelta = delta_var/(double)(X_POINTS-1);
    double ydelta = (max_eval - min_eval)/(double)(Y_POINTS-1);

    for(long long ict = 0; ict < X_POINTS; ict ++){
      input_var = min_eval_range + (double)(ict)*gdelta;
      _function.execute();
      for(long long jct = 0; jct < Y_POINTS; jct ++){
        if((output_var >= (min_eval-ydelta/2.0 + ydelta*(double)jct))
            && (output_var <=  (min_eval + ydelta/2.0 + ydelta*(double)jct))){
          gchars[ict][jct] = '*';
        }
        else{
          gchars[ict][jct] = ' ';
        }
      }
    }

    for(long long jct = Y_POINTS-1; jct >= 0; jct --){
      input_var = min_eval_range + (double)(jct)*gdelta;
      _function.execute();
      char char_array[10];
      char_array[0] = '\0';
      double var = min_eval+ydelta*(double)jct;
      sprintf(char_array,"%6.2e",var);

      s <<  char_array  << "\t|";
      for(long long ict = 0; ict < X_POINTS; ict ++){
        s << gchars[ict][jct];
      }
      s << std::endl;
    }s << "\t\t--------------------------------------------------------------" << std::endl;
    s << "\t\tPlot of function across mesh range." << std::endl << std::endl << std::endl;

  }


  /*****************************************************************************/
  Element_Density_Function::Element_Density_Function(
      const std::string & funcBody,
      const std::string & app_dir,
      std::stringstream & error_string
      ):
    _funcBody(funcBody),
    application_direction(app_dir),
    _function(2), //passing in by reference coord(input value) and field(output value)
    integrated(false),
    integral_total(0.0),
    max_eval(0.0),
    min_eval(0.0),
    min_eval_range(0.0),
    max_eval_range(0.0)
  {
    _function.addVar("double", "coord");
    _function.addVar("double", "field");
    bool success = _function.addBody(_funcBody);
    if (!success) {
      error_string << "User_Defined_Element_Density_Function::User_Defined_Element_Density_Function: "
        << "function body error: " << _function.getErrors();
    }
    running_sum = NULL;
  }


  /*****************************************************************************/
  Element_Density_Function::~Element_Density_Function()
    /*****************************************************************************/
  {
    if(running_sum) delete [] running_sum;
    running_sum = NULL;
  }

  /*****************************************************************************/
  void Element_Density_Function::Integrate(double start_var, double end_var,std::stringstream & es)
    /*****************************************************************************/
  {
    min_eval_range = start_var;
    max_eval_range = end_var;

    long long num_points = running_sum_length + 1;
    running_sum = new double[num_points];
    double input_var = 0.;
    double * ivp = & input_var;
    _function.varAddrFill(0, ivp);
    double output_var = 0;
    double * ovp = & output_var;
    _function.varAddrFill(1, ovp);
    integral_total = 0;
    running_sum[0] = integral_total;
    double delta_var = max_eval_range-min_eval_range;
    double delta = delta_var/running_sum_length;

    //initialize eval numbers;

    input_var = min_eval_range;
    _function.execute();
    min_eval = max_eval = output_var;

    // cumulative trapezoidal rule integration
    for(long long i = 1; i < num_points; i ++){

      input_var = min_eval_range + (double)(i-1)*delta;
      _function.execute();
      if(output_var <=0.){
        es << "Element_Density_Function::Integrate(): "
          << "User defined run-time compled function evaluates to 0. or negative value. \n"
          << "User defined run time functions used for mesh grading must evaluate\n"
          << "to a positive value across the range of the mesh.";
      }
      if(output_var > max_eval)max_eval = output_var;
      if(output_var < min_eval)min_eval = output_var;
      integral_total += output_var;

      input_var = min_eval_range + (double)(i)*delta;
      _function.execute();
      if(output_var <=0.)
        es << "Element_Density_Function::Integrate(): "
          << "User defined run-time compled function evaluates to 0. or negative value. \n"
          << "User defined run time functions used for mesh grading must evaluate\n"
          << "to a positive value across the range of the mesh.";
      if(output_var > max_eval)max_eval = output_var;
      if(output_var < min_eval)min_eval = output_var;
      integral_total += output_var;
      running_sum[i] = integral_total/2.0;
    }
    integral_total = integral_total/2.0;
    for(long long i = 0; i < num_points; i ++){
      running_sum[i]/=integral_total;
    }
  }


  /*****************************************************************************/
  double Element_Density_Function::Interpolate(double incoming_var, std::stringstream & es)
    /*****************************************************************************/
  {
    long long num_points = running_sum_length + 1;
    for(long long i = 0; i < num_points; i ++){
      if(running_sum[i] == incoming_var)return (double)i/(double)running_sum_length;
      if(running_sum[i+1] == incoming_var)return (double)(i+1)/(double)running_sum_length;
      if(incoming_var > running_sum[i] && incoming_var < running_sum[i+1]){
        //interpolate
        double imin = (double)running_sum[i];double imax = (double)running_sum[i+1];
        double delta = incoming_var-imin;
        double full_delta = imax-imin;
        return((double)i/(double)running_sum_length + (delta/full_delta)/(double)running_sum_length);
      }
    }
    es <<"User_Defined_Element_Density::Interpolate: "
      << "interpolant not found" ;
    return (-1.);
  }


}// end namespace
