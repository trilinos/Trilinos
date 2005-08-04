//
// @HEADER
// ***********************************************************************
// 
//                           Capo Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Capo_PARAMETER_LIST_H
#define Capo_PARAMETER_LIST_H

/************************************************************ 
File:      Capo_Parameter_List.hpp
Purpose:   The Parameter List object contains all the necessary
           parameters for specific Capo algorithms.
Date:      6-10-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace CAPO {
  
  class Integrator;
  
  /** \brief A list of Parameters to be set to influence the
      Capo algorithms.
  */
  class Parameter_List {
  public:
    //! Constructor
    Parameter_List();

    //! Destructor
    ~Parameter_List();

    //! Returns the maximum allowed number of continuation steps.
    int get_MaxOuterIts();
    //! Returns the number of inner iterations allowed.
    int get_MaxInnerIts();
    //! Returns the corrector tolerance.
    double get_tol();
    //! Returns the printing output level.
    int get_printproc();
    //! Returns the initial lambda stepsize.
    double get_lambda_stepsize();
    //! Returns the maximum allowed lambda step size.
    double get_lambda_max();
    //! Returns the minimum allowed lambda step size.
    double get_lambda_min();
    //! Returns the number of inner iterations which signals
    //! an allowed increase in the lambda step size.
    int get_lambda_extend_tol();
    //! Returns true if the problem exhibits a periodic solution.
    bool Periodic();
    //! Returns true if arc-length continuation is turned on.
    bool Arc_Length();
    //! Returns the number of subspace iterations to be performed.
    int get_SubspaceIterations();
    //! Returns the floquet multiplier tolerance.
    double get_FloquetTolerence();
    //! Returns the number off additional vectors in the "unstable" subspace.
    int get_NumberXtraVecsSubspace();

    //! Set integer parameters.
    void set_param(string& name, int value);
    //! Set boolean parameters.
    void set_param(string& name, bool value);
    //! Set double parameters.
    void set_param(string& name, double value);


  private:
    //! Number of continuation steps to take.  Default = 20.
    int MaxOuterIts;

    //! Number of steps for the corrector.  Default = 10.
    int MaxInnerIts;

    //! Corrector Solve Tolerance. Default = 1.0e-5.
    double tol;

    //! Sets the printing output level. Default = 1.
    int printproc;

    //! Continuation parameter stepsize. Default = 0.1.
    double lambda_stepsize;

    //! Tolerance for extending the lambda_stepsize. Default = 3.
    int lambda_extend_tol;

    //! Maximum lambda stepsize. Default = .5.
    double lambda_max;

    //! Minimum lambda stepsize. Default = .00001.
    double lambda_min;

    //! Searching for a Periodic Solution if true
    //! otherwise looking for a fixed point.  Default=false.
    bool EnablePeriodicity;

    //! Turn Arclength continuation on or off.  Default = fales.
    bool EnableArclength;

    //! Number of subspace iterations.  Default = 10.
    int SubspaceIterations;

    //! Cut off for a floquet multiplier to be determined large enough
    //! to be in the "unstable" realm. Default = .5.
    double FloquetTolerence;

    //! Allowed number of extra vectors in the unstable subspace to 
    //! help speed up convergence. Default = 4.
    int NumberXtraVecsSubspace;


  };

}

#endif
