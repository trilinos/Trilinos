// $Id$ 
// $Source$ 

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

#ifndef NOX_TSF_SHARED_OPERATOR_H
#define NOX_TSF_SHARED_OPERATOR_H

#include "NOX_Common.H"
#include "TSFLinearOperator.hpp"

//class TSFExtended::LinearOperator<double>;

namespace NOX {
namespace TSF {

class Group;			// forward declaration

  /*! \brief Shared operator object for NOX::TSF::Group objects.
   * 
   * Due to the size of certain operators (i.e. the Jacobian and 
   * Preconditioning Matrices), we cannot afford to have
   * multiple copies. Instead we implement a shared operator class that
   * that all groups use. 
   *
   * This class is not needed, since TSFExtended::LinearOperator already
   * shares memory. KL 13 July 2004. 
   */

class SharedOperator {

public:

  //! Constructor
  SharedOperator(TSFExtended::LinearOperator<double> A);

  //! Destructor
  ~SharedOperator();

  //! Get a non-const reference to the underlying jacobian matrix (a TSFExtended::LinearOperator)
  TSFExtended::LinearOperator<double> getLinearOperator(const Group* newowner);
 
  //! Return a const reference to the Jacobian
  const TSFExtended::LinearOperator<double> getLinearOperator() const;
  
  //! Return true if grp is the owner of the Jacobian
  bool isOwner(const Group* grp) const;


protected:



  //! TSFExtended::LinearOperator 
  TSFExtended::LinearOperator<double> linearOp;


  //! Pointer to const group that owns Jacobian
  const Group* owner;
};
} // namespace TSF
} // namespace NOX

#endif
