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


/************************************************************ 
File:      Capo_Parameter_List.cpp
Purpose:   The Parameter List object contains all the necessary
           parameters for specific Capo algorithms.
Date:      6-10-05
Author:    Joseph Simonis
**************************************************************/

/**** Includes ****/
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Capo_Parameter_List.hpp"

using namespace CAPO;

//-----------------------------------------------------------------
// Function      : Parameter_List::Parameter_List
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
Parameter_List::Parameter_List()
{
  // Set the parameters to defaults:
  MaxOuterIts = 20;
  MaxInnerIts = 10;
  tol = 1.0e-5;
  printproc = 1;
  lambda_stepsize = .1;
  lambda_extend_tol = 3;
  lambda_max = .5;
  lambda_min = .00001;
  EnablePeriodicity = false;
  EnableArclength = false;

  // Parameters for Newton-Picard Gauss-Seidel
  SubspaceIterations = 10;
  FloquetTolerence = .5;
  NumberXtraVecsSubspace = 4;
  
  // Parameters for Recursive Projection Method
  Mplus2tol = .1;
  ModifiedNewtFlag = 0;
  UpdateBasisFreq =2;
  cerr << "Finished with Parameter_List constructor." << endl;
}

//-----------------------------------------------------------------
// Function      : Parameter_List::~Parameter_List
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
Parameter_List::~Parameter_List()
{
}

//-----------------------------------------------------------------
// Function      : Parameter_List::set_param
// Purpose       : method for setting parameters of type int
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Parameter_List::set_param(string& name, int value)
{
  bool found = false;
  if (name=="SubspaceIterations")
    {
      SubspaceIterations = value;
      found = true;
    }
  if (name=="NumberXtraVecsSubspace")
    {
      NumberXtraVecsSubspace = value;
      found = true;
    }
  if (name=="ModifiedNewtFlag")
    {
      ModifiedNewtFlag = value;
      found = true;
    }
  if (name=="UpdateBasisFreq")
    {
      UpdateBasisFreq = value;
      found = true;
    }
  if (name=="MaxOuterIts")
    {
      MaxOuterIts = value;
      found = true;
    }

  if (name=="MaxInnerIts")
    {
      MaxInnerIts = value;
      found = true;
    }
  if (name=="printproc")
    {
      printproc = value;
      found = true;
    }
  if (name=="lambda_extend_tol")
    {
      lambda_extend_tol = value;
      found = true;
    }

  if (!found)
    cout << "Parameter " << name << " not found!" << endl;
}

//-----------------------------------------------------------------
// Function      : Parameter_List::set_param
// Purpose       : method for setting parameters of type bool
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Parameter_List::set_param(string& name, bool value)
{
  bool found = false;
  if (name=="EnableArclength")
    {
      EnableArclength = value;
      found = true;
    }
  if (name=="EnablePeriodicity")
    {
      EnablePeriodicity = value;
      found = true;
    }
  if (!found)
    cout << "Parameter " << name << " not found!" << endl;


}

//-----------------------------------------------------------------
// Function      : Parameter_List::set_param
// Purpose       : method for setting parameters of type double
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
void Parameter_List::set_param(string& name, double value)
{
  bool found = false;
  if (name=="FloquetTolerence")
    {
      FloquetTolerence = value;
      found = true;
    }
  if (name=="Mplus2tol")
    {
      Mplus2tol = value;
      found = true;
    }
  if (name=="tol")
    {
      tol = value;
      found = true;
    }
  if (name=="lambda_stepsize")
    {
      lambda_stepsize = value;
      found = true;
    }
  if (name=="lambda_max")
    {
      lambda_max = value;
      found = true;
    }
  if (name=="lambda_min")
    {
      lambda_min = value;
      found = true;
    }
  if (!found)
    cout << "Parameter " << name << " not found!" << endl;
}

//-----------------------------------------------------------------
// Function      : Parameter_List::get_"param"
// Purpose       : the following functions return the values
//                 of specified parameters.
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 06/10/05
//------------------------------------------------------------------
int Parameter_List::get_MaxOuterIts()
{return MaxOuterIts;}
int Parameter_List::get_MaxInnerIts()
{return MaxInnerIts;}
double Parameter_List::get_tol()
{return tol;}
int Parameter_List::get_printproc()
{return printproc;}
double Parameter_List::get_lambda_stepsize()
{return lambda_stepsize;}
double Parameter_List::get_lambda_max()
{return lambda_max;}
double Parameter_List::get_lambda_min()
{return lambda_min;}

int Parameter_List::get_SubspaceIterations()
{ return SubspaceIterations;}
double Parameter_List::get_FloquetTolerence()
{ return FloquetTolerence;}
int Parameter_List::get_NumberXtraVecsSubspace()
{ return NumberXtraVecsSubspace;}
int Parameter_List::get_lambda_extend_tol()
{return lambda_extend_tol;}

double Parameter_List::get_Mplus2tol()
{ return Mplus2tol;}
int Parameter_List::get_ModifiedNewtFlag()
{ return ModifiedNewtFlag;}
int Parameter_List::get_UpdateBasisFreq()
{ return UpdateBasisFreq;}

//-----------------------------------------------------------------
// Function      : Parameter_List::Periodic
// Purpose       : Returns the value of EnablePeriodicty
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/01/05
//------------------------------------------------------------------
bool Parameter_List::Periodic()
{ return EnablePeriodicity;}


//-----------------------------------------------------------------
// Function      : Parameter_List::Arc_Length
// Purpose       : Returns the value of EnableArclength
// Special Notes :
// Scope         : public
// Creator       : J. Simonis, SNL
// Creation Date : 07/01/05
//------------------------------------------------------------------
bool Parameter_List::Arc_Length()
{ return EnableArclength;}
