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

#ifndef STOKHOS_MP_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_MP_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_MPPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! Factory for generating stochastic Galerkin preconditioners
  class MPPreconditionerFactory {
  public:

    //! Constructor
    MPPreconditionerFactory(
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~MPPreconditionerFactory() {}

    //! Build preconditioner operator
    virtual Teuchos::RCP<Stokhos::MPPreconditioner> 
    build(
      const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      int num_mp_blocks,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& mp_map);

  protected:

    //! Build preconditioner factory for each point
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> 
    buildPointPreconditionerFactory();

  private:
    
    //! Private to prohibit copying
    MPPreconditionerFactory(const MPPreconditionerFactory&);
    
    //! Private to prohibit copying
    MPPreconditionerFactory& operator=(const MPPreconditionerFactory&);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

  }; // class MPPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_MP_PRECONDITIONER_FACTORY_HPP
