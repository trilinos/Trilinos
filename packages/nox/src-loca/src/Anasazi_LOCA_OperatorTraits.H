// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef ANASAZI_LOCA_OPERATORTRAITS_H
#define ANASAZI_LOCA_OPERATORTRAITS_H

// Anasazi includes
#include "AnasaziOperatorTraits.hpp"

// NOX/LOCA includes
#include "NOX_Abstract_MultiVector.H"
#include "LOCA_AnasaziOperator_AbstractStrategy.H"

namespace Anasazi {

  /*!
   * \brief Implementation of the Anasazi::OperatorTraits for LOCA operators.
   */
  /*!
   * This class provides an implementation of the Anasazi::OperatorTraits class
   * using LOCA operators.  The operator is represented by the 
   * LOCA::AnasaziOperator::AbstractStrategy class which serves as a base class
   * for many different types of operators, including shift-invert and
   * Generalized Cayley transformations.
   */
  template <>
  class OperatorTraits<double,
		       NOX::Abstract::MultiVector,
		       LOCA::AnasaziOperator::AbstractStrategy>
  {
  public:

    //! Synony for double
    typedef double ScalarType;

    //! Synonym for NOX::Abstract::MultiVector
    typedef NOX::Abstract::MultiVector MV;

    //! Synonym for LOCA::AnasaziOperator::AbstractStrategy
    typedef LOCA::AnasaziOperator::AbstractStrategy OP;
    
    //@{ Matrix/Operator application method.
    
    //! Application method which performs operation <b>y = Op*x</b>
    /*!
     *\return Status of the operation.
     */
    static void Apply ( const OP& Op, const MV& x, MV& y )
    { 
      Op.apply(x,y);
    };
    
    //@}
    
  }; // class MultiVecTraits
  
} // namespace Anasazi

#endif 
