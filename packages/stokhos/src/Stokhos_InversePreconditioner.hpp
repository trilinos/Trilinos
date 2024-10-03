// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_INVERSEPRECONDITIONER_HPP
#define STOKHOS_INVERSEPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class InversePreconditioner : 
    public Stokhos::Operator<ordinal_type,double> {
  public:

    //! Constructor 
    InversePreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type,double> & A_) : A(A_) {}
    
    //! Destructor
    virtual ~InversePreconditioner() {}
    
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type m) const {
      Teuchos::RCP< Teuchos::SerialDenseMatrix<ordinal_type, value_type> > AA, UU, RR;
      AA = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (Teuchos::Copy,A));
      UU = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (Teuchos::Copy,Result));
      RR = Teuchos::rcp(new Teuchos::SerialDenseMatrix<ordinal_type, value_type> (Teuchos::Copy,Input));

      // Setup solver
      Teuchos::SerialDenseSolver<ordinal_type, value_type> solver;
      solver.setMatrix(AA);
      solver.setVectors(UU, RR);
      //Solve A*Result=Input
      if (solver.shouldEquilibrate()) {
         solver.factorWithEquilibration(true);
         solver.equilibrateMatrix();
      }
      solver.solve();
      
      for (ordinal_type i=0; i<A.numRows(); i++)
	Result(i,0)=(*UU)(i,0);
      
      return 0;
    }
   
  protected:
    const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & A;
  }; // class InversePreconditioner

} // namespace Stokhos

#endif // STOKHOS_INVERSEPRECONDITIONER_HPP

