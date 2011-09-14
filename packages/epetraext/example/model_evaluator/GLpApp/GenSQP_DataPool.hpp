/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef GENSQP_DATAPOOL_H
#define GENSQP_DATAPOOL_H

#include "GenSQP_Vector.hpp"

/** \class GenSQP::DataPool
    \brief Provides the interface to a generic data pool.

    The only member function is computeAll, responsible for the computation of all data that can be
    stored and reused within one (or more) SQP iteration. The data is entirely problem specific.
*/


namespace GenSQP {

class DataPool {
public:

  virtual ~DataPool() {}
  
  /** \brief Recompute all stored quantities.
      \param x [in]  - Current SQP iterate vector.
    
      \return None. 
    
      \par Detailed Description:

      Interface function that evaluates and stores problem-specific quantities
      that can be reused within one (or more) SQP iteration.
  
      \note The GenSQP::Vector input parameters can be recast into user-accessible
      quantities by using the following syntax:\n
      <tt>
      Teuchos::RefCountPtr<const umDVec> ex =
        (Teuchos::dyn_cast<GenSQP::SledgeVector>(const_cast<GenSQP::Vector&>(x))).getVector();
      </tt>
  */
  virtual void computeAll( const Vector &x ) = 0;

}; // class DataPool

} // namespace GenSQP

#endif
