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

#ifndef LOCA_SAVEEIGENDATA_DEFAULTSTRATEGY_H
#define LOCA_SAVEEIGENDATA_DEFAULTSTRATEGY_H

#include "LOCA_SaveEigenData_AbstractStrategy.H" // base class

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

  namespace SaveEigenData {

    //! Default strategy for saving eigenvector/value data
    /*!
     * This class implements a default strategy for saving eigenvectors
     * and eigenvalues that does nothing and exists so the LOCA::Stepper
     * always has an object to pass eigen data to.
     */
    class DefaultStrategy : public LOCA::SaveEigenData::AbstractStrategy {

    public:

      //! Constructor
      DefaultStrategy(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);
      
      //! Destructor
      virtual ~DefaultStrategy();

      //! Save eigenvalues/eigenvectors
      /*!
       * The implementation here does nothing and always returns
       * NOX::Abstract::Group::Ok.
       */
      virtual NOX::Abstract::Group::ReturnType 
      save(Teuchos::RCP< std::vector<double> >& evals_r,
	   Teuchos::RCP< std::vector<double> >& evals_i,
	   Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
	   Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i);

    private:

      //! Private to prohibit copying
      DefaultStrategy(const DefaultStrategy&);

      //! Private to prohibit copying
      DefaultStrategy& operator = (const DefaultStrategy&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class DefaultStrategy

  } // Namespace SaveEigenData

} // Namespace LOCA

#endif
