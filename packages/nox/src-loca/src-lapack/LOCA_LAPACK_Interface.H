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

#ifndef LOCA_LAPACK_INTERFACE_H
#define LOCA_LAPACK_INTERFACE_H

#include "NOX_LAPACK_Vector.H"
#include "NOX_LAPACK_Interface.H"

// Forward declarations
namespace LOCA {
  class ParameterVector;
}

namespace LOCA {
  namespace LAPACK {

    /*! 
      \brief Virtual interface for LOCA::LAPACK::Group to provide a
      link to the external code for RHS and Jacobian fills. 
    */
    /*!
      This interface is derived from the NOX::LAPACK::Interface and 
      and additionally provides a method for setting problem parameters.
    */
    class Interface : public NOX::LAPACK::Interface {

    public:
 
      //! Constructor
      Interface() {}

      //! Destructor
      virtual ~Interface() {}

      //! Set the problem parameters
      virtual void setParams(const ParameterVector& p) = 0;

      //! Call user's own print routine for vector-parameter pair
      virtual void printSolution(const NOX::LAPACK::Vector& x_,
                                 const double conParam) {}

      //! Compute mass matrix.  Returns true if computation was successful.
      virtual bool computeMass(NOX::LAPACK::Matrix& M, 
			       const NOX::LAPACK::Vector & x) { return false; }

      //! Projects solution to a few scalars for multiparameter continuation
      /*! Default implementation is the max norm 
       */
      virtual void projectToDraw(const NOX::LAPACK::Vector& x,
				 double *px) const {
	px[0] = x.norm(NOX::Abstract::Vector::MaxNorm);
      }

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const { return 1; }

    };
  } // namespace LAPACK
} // namespace LOCA

#endif
