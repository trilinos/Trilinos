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

#ifndef LOCA_EIGENSOLVER_DGGEVSTRATEGY_H
#define LOCA_EIGENSOLVER_DGGEVSTRATEGY_H

#include "LOCA_Eigensolver_AbstractStrategy.H" // base class

// Forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Eigensolver {

    //! LAPACK DGGEV eigensolver strategy
    /*!
     * This class implements an eigensolver strategy using LAPACK routines
     * DGGEV and DGEEV
     */
    class DGGEVStrategy : public LOCA::Eigensolver::AbstractStrategy {

    public:

      //! Constructor
      DGGEVStrategy(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);
      
      //! Destructor
      virtual ~DGGEVStrategy();

      //! Compute eigenvalues/eigenvectors
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(
	      NOX::Abstract::Group& group,
	      Teuchos::RCP< std::vector<double> >& evals_r,
	      Teuchos::RCP< std::vector<double> >& evals_i,
	      Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
	      Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      DGGEVStrategy(const DGGEVStrategy&);

      //! Private to prohibit copying
      DGGEVStrategy& operator = (const DGGEVStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parsed top-level parameter list
      Teuchos::RCP<LOCA::Parameter::SublistParser> topParams;

      //! Eigensolver parameter list
      Teuchos::RCP<Teuchos::ParameterList> eigenParams;

      //! The number of requested eigenvalues
      int nev;

      //! Which eigenvalues are of interest
      string which;

    }; // Class DGGEVStrategy

  } // Namespace Eigensolver

} // Namespace LOCA

#endif
