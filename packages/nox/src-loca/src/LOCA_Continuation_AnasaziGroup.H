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

#ifndef LOCA_CONTINUATION_ANASAZIGROUP_H
#define LOCA_CONTINUATION_ANASAZIGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

namespace LOCA { 

  namespace Continuation { 

    /*! 
     * \brief %LOCA abstract interface for continuation, derived from the 
     * NOX::Abstract::Group.  This abstract class provides 
     * the interface necessary to perform continuation, i.e., compute 
     * families of solutions to \f$ F(x,p) = 0 \f$.
    */
    /*! 
     * Concrete implemenations of this interface must provide implementations 
     * of all of the methods in the NOX::Abstract::Group interface as well as 
     * the additional interface defined here.
    */

    class AnasaziGroup : public virtual LOCA::Continuation::AbstractGroup {

    public:
  
      //! Default constructor.
      AnasaziGroup();

      //! Copy constructor.
      AnasaziGroup(const AnasaziGroup& source, 
		   NOX::CopyType type = NOX::DeepCopy);
 
      //! Destructor.
      virtual ~AnasaziGroup();

      //! Assignment operator
      virtual AnasaziGroup& operator=(const AnasaziGroup& source);

      //! Computes leading eigenvalues/eigenvectors using Anasazi package
      /*!
       * To use Anasazi for computing eigenvalues, NOX/LOCA must be
       * configured with --enable-loca-anasazi
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(NOX::Parameter::List& params);

     NOX::Abstract::Group::ReturnType  
     applyAnasaziOperator (NOX::Parameter::List &params, 
                           const NOX::Abstract::Vector &input, 
                           NOX::Abstract::Vector &result) const;

    protected:

      /*! 
       * \brief Step counter for option of computing eigenvalues every n 
       * continuation steps
       */
      int eigenvalCounter;
      bool hasMassMatrix;

    };
  } // namespace Continuation
} // namespace LOCA

#endif
