// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_ANASAZIGROUP_H
#define LOCA_CONTINUATION_ANASAZIGROUP_H

#include "LOCA_Continuation_AbstractGroup.H"    // base class

// forward declarations
#ifdef HAVE_LOCA_ANASAZI
namespace Anasazi {
  class LOCAMultiVec;
}
#endif

namespace LOCA { 

  namespace Continuation { 

    /*! 
     * \brief %LOCA eigenvalue driver using Anasazi
     */
    /*! 
     * This class provides a concrete implementation of the 
     * computeEigenvalues() method of the LOCA::Continuation::AbstractGroup
     * using Anasazi, the generic Trilinos eigensolver.  Eigenvalues of the
     * Jacobian and generalized eigenvalues of the Jacobian/mass matrix
     * are supported using shift-invert, as well as generalized Cayley
     * transformations.
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
      //! Driver for eigensolver, including preprocessing of parameters,
      //! generating initial vector, calling eigensolver, and post-processing
      //! of the results.
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(NOX::Parameter::List& params);

//       //! Routine called back to by Anasazi to perform one application
//       //! of the operator. Contains logic for shift-and-invert and
//       //! the generalized Cayley transformation.
//       NOX::Abstract::Group::ReturnType  
//       applyAnasaziOperator (NOX::Parameter::List &params, 
//                             const NOX::Abstract::Vector &input, 
//                             NOX::Abstract::Vector &result) const;

#ifdef HAVE_LOCA_ANASAZI
      //! Save eigenvectors and eigenvalues for later processing
      /*!
       * A derived group is free to define this method as 
       * appropriate.  Typically they are saved by the application code
       * in an output file for later processing/restarting.
       *
       * The default implementation here does nothing.
       */
      virtual void
      saveEigenVectors(int nev, const std::vector<double>& evals, 
		       const Anasazi::LOCAMultiVec& evecs ) const;
#endif

    protected:

      /*! 
       * \brief Step counter for option of computing eigenvalues every n 
       * continuation steps
       */
      int eigenvalCounter;

//       //! Flag whether the system has a mass matrix, otherwise the identity matrix is used.
//       bool hasMassMatrix;

//       //! Flag whether the Cayley transformation is used, otherwise Shift-and-Invert is used.
//       bool cayley;

//       //! The shift (pole) parameter in the Cayley transformation.
//       double sigma;

//       //! The 'zero' parameter in the Cayley transformation.
//       double mu;

//       //! The shift (pole) parameter in the shift-and-Invert transformation.
//       double shift_;
    };
  } // namespace Continuation
} // namespace LOCA

#endif
