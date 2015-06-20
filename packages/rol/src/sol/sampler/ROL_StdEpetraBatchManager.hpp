// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_STDEPETRABATCHMANAGER_HPP
#define ROL_STDEPETRABATCHMANAGER_HPP

#include "ROL_EpetraBatchManager.hpp"
#include "ROL_StdVector.hpp"

namespace ROL {

template<class Real> 
class StdEpetraBatchManager : public EpetraBatchManager<Real> {
public:
  StdEpetraBatchManager(Teuchos::RCP<Epetra_Comm> &comm) : EpetraBatchManager<Real>(comm) {}
  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    Teuchos::RCP<std::vector<Real> > input_ptr = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<StdVector<Real> >(input)).getVector());
    Teuchos::RCP<std::vector<Real> > output_ptr = Teuchos::rcp_const_cast<std::vector<Real> >(
      (Teuchos::dyn_cast<StdVector<Real> >(output)).getVector());
    int dim_i = input_ptr->size();
    int dim_o = output_ptr->size();
    if ( dim_i != dim_o ) {
      std::cout << "StdEpetraBatchManager: DIMENSION MISMATCH ON RANK " 
                << EpetraBatchManager<Real>::batchID() << "\n";
    }
    else {
      EpetraBatchManager<Real>::sumAll(&(*input_ptr)[0],&(*output_ptr)[0],dim_i);
    }
  }
};

}

#endif
