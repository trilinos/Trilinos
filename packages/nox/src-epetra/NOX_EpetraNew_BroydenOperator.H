//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef _NOX_EPETRANEW_BROYDENOPERATOR_H
#define _NOX_EPETRANEW_BROYDENOPERATOR_H

#include "NOX.H"
#include "NOX_Common.H"
#include "NOX_Epetra.H"
#include "NOX_Parameter_List.H"
#include "NOX_Parameter_PrePostOperator.H"

#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

namespace NOX {

namespace EpetraNew {

/*! \brief Pure virtual class for creating an Operator for NOX.
  This class derives from the Epetra_Operator class with a single
  additional method: compute().    This allows the user to supply their own
  preconditioner or jacobian operators.
 */
class BroydenOperator : public NOX::Parameter::PrePostOperator,
                        public NOX::EpetraNew::Interface::Preconditioner {

 public:

  //! Constructor
  BroydenOperator(NOX::Parameter::List & nlParams, Epetra_Vector & solnVec,
                  Epetra_CrsMatrix & mat);

  //! Copy Constructor
  BroydenOperator(const BroydenOperator &);

  //! Pure virtual destructor
  virtual ~BroydenOperator();

  Epetra_CrsMatrix & getBroydenMatrix()
  { return crsMatrix; };

  // 
  virtual bool computePreconditioner( const Epetra_Vector &, 
                                      NOX::Parameter::List * params = 0 );

  //
  virtual void runPostIterate( const NOX::Solver::Generic & solver);

  virtual NOX::Parameter::Arbitrary * clone() const;

  virtual const string & getType() const; 

 protected:

  NOX::Epetra::Vector * updateVectorPtr;
  NOX::Epetra::Vector & updateVector;

  NOX::Epetra::Vector * broydenVecPtr;
  NOX::Epetra::Vector * residualVecPtr;

  Epetra_CrsMatrix & crsMatrix;

  std::string myType;
};
} // namespace EpetraNew
} // namespace NOX

#endif /* _NOX_EPETRANEW_BROYDENOPERATOR */
