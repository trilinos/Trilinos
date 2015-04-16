/*
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
*/


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

