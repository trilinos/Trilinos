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

#ifndef LOCA_THYRA_SAVEDATASTRATEGY_H
#define LOCA_THYRA_SAVEDATASTRATEGY_H

#include "Teuchos_RCP.hpp" 
#include "NOX_Abstract_Vector.H" 
#include "LOCA_Abstract_Iterator.H"

namespace LOCA {

  namespace Thyra {
    
    //! Abstract interface class strategies to save continuation data
    class SaveDataStrategy {

    public:

      //! Constructor
      SaveDataStrategy() {}

      //! Destructor
      virtual ~SaveDataStrategy() {}

      //! Save solution
      /*!
       * \param x [in] Solution vector
       * \param p [in] Parameter value
       * \returns ReturnType code indicating success or failure
       */
      virtual void
      saveSolution(const NOX::Abstract::Vector& x, double p) {}

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The default implementation to empty.
       */
      virtual void 
      preProcessContinuationStep(
			    LOCA::Abstract::Iterator::StepStatus stepStatus) {}

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.  The default implementation to empty.
       */
      virtual void 
      postProcessContinuationStep(
			    LOCA::Abstract::Iterator::StepStatus stepStatus) {}

      //! Projects solution to a few scalars for multiparameter continuation
      /*! 
       * This method is called every time a solution is saved by the 
       * multiparameter continuation code MF for later visualization
       * and should project the solution vector down to a few scalars.
       * The array \c px will be preallocated to the proper length
       * given by projectToDrawDimension().
       */
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
				 double *px) const {}

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const { return 0; }

    private:

      //! Private to prohibit copying
      SaveDataStrategy(const SaveDataStrategy&);

      //! Private to prohibit copying
      SaveDataStrategy& operator = (const SaveDataStrategy&);

    }; // Class SaveDataStrategy

  } // Namespace Thyra

} // Namespace LOCA

#endif // LOCA_THYRA_SAVEDATASTRATEGY_H
