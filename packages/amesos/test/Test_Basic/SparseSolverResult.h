/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
*/

#ifndef SPARSESOLVERRESULT
#define SPARSESOLVERRESULT
#include "Epetra_Object.h"
const double UnUsedDbl = 13e130;

class SparseSolverResult : public Epetra_Object { 
  
 public:
  SparseSolverResult() :  
    first_time( UnUsedDbl ), middle_time( UnUsedDbl ), 
    last_time( UnUsedDbl ), total_time( UnUsedDbl ),
    error( UnUsedDbl), residual(UnUsedDbl), 
    Anorm( UnUsedDbl ), Xnorm( UnUsedDbl ), Bnorm( UnUsedDbl ) 
  { 
  ; } ; 
  ~SparseSolverResult(){};

  void Set_First_Time ( double First_Time_in ) { first_time = First_Time_in; } ; 
  inline double Get_First_Time ( ) { return first_time ; } ; 

  void Set_Middle_Time ( double Middle_Time_in ) { middle_time = Middle_Time_in; } ; 
  inline double Get_Middle_Time ( ) { return middle_time ; } ; 

  void Set_Last_Time ( double Last_Time_in ) { last_time = Last_Time_in; } ; 
  inline double Get_Last_Time ( ) { return last_time ; } ; 

  void Set_Total_Time ( double Total_Time_in ) { total_time = Total_Time_in; } ; 
  inline double Get_Total_Time ( ) { return total_time ; } ; 

  void Set_Bnorm ( double bnorm_in ) { Bnorm = bnorm_in; } ; 
  inline double Get_Bnorm ( ) { return Bnorm ; } ; 
  void Set_Xnorm ( double xnorm_in ) { Xnorm = xnorm_in; } ; 
  inline double Get_Xnorm ( ) { return Xnorm ; } ; 
  void Set_Residual ( double residual_in ) { residual = residual_in; } ; 
  inline double Get_Residual ( ) { return residual ; } ; 
  void Set_Error ( double error_in ) { error = error_in; } ; 
  inline double Get_Error ( ) { return error; } ; 
  void Set_Anorm ( double anorm_in ) { Anorm = anorm_in; } ; 
  inline double Get_Anorm ( ) { return Anorm; } ; 

  virtual void Print(std::ostream & os) const;
  virtual void PrintSummary(std::ostream & os) const;


 private:
  double first_time ;
  double middle_time ;
  double last_time ;
  double total_time ;
  double error ;
  double residual ;
  double Anorm ;
  double Xnorm ;
  double Bnorm ;
} ;

#endif
