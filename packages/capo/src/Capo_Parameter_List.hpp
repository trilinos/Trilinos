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
  
  /**** Parameter_List Class ****/
  class Parameter_List {
  public:
    //! Constructor
    Parameter_List();

    //! Destructor
    ~Parameter_List();

    //! Member functions
    int get_MaxOuterIts();
    int get_MaxInnerIts();
    double get_tol();
    int get_printproc();
    double get_lambda_stepsize();
    bool Periodic();
    bool Arc_Length();

    //! Member functions that will be used by Newton-Picard Gauss-Seidel
    int get_SubspaceIterations();
    double get_FloquetTolerence();
    int get_NumberXtraVecsSubspace();

    //! Member functions that will be used by Recursive Projection Method
    double get_Mplus2tol();
    int get_ModifiedNewtFlag();
    int get_UpdateBasisFreq();

    //! Functions to set parameters
    void set_param(string& name, int value);
    void set_param(string& name, bool value);
    void set_param(string& name, double value);


  private:

    //! Shared Parameters
    int MaxOuterIts;
    int MaxInnerIts;
    double tol;
    int printproc;
    double lambda_stepsize;
    bool EnablePeriodicity;
    bool EnableArclength;


    //! Parameters for Newton-Picard Gauss-Seidel
    int SubspaceIterations;
    double FloquetTolerence;
    int NumberXtraVecsSubspace;

    //! Parameters for Recursive Projection Method
    double Mplus2tol;
    int ModifiedNewtFlag;
    int UpdateBasisFreq;

  };

}

#endif
