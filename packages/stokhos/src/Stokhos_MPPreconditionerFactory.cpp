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

#include "Stokhos_MPPreconditionerFactory.hpp"
#include "Stokhos_MPBlockDiagonalPreconditioner.hpp"
#include "Stokhos_MPMeanBasedPreconditioner.hpp"
#include "Stokhos_PreconditionerFactory.hpp"
#include "Teuchos_Assert.hpp"

Stokhos::MPPreconditionerFactory::
MPPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  params(params_)
{
}

Teuchos::RCP<Stokhos::MPPreconditioner> 
Stokhos::MPPreconditionerFactory::
build(const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      int num_mp_blocks,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& mp_map)
{
  Teuchos::RCP<Stokhos::MPPreconditioner> mp_prec;
  std::string prec_method = params->get("Preconditioner Method", 
					"Block Diagonal");
  if (prec_method == "Block Diagonal") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildPointPreconditionerFactory();
    mp_prec = Teuchos::rcp(new Stokhos::MPBlockDiagonalPreconditioner(
			     mp_comm, num_mp_blocks, 
			     base_map, mp_map, prec_factory, 
			     params));
  }
  else if (prec_method == "Mean-based") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildPointPreconditionerFactory();
    mp_prec = Teuchos::rcp(new Stokhos::MPMeanBasedPreconditioner(
			     mp_comm, num_mp_blocks, 
			     base_map, mp_map, prec_factory, 
			     params));
  }
  else if (prec_method == "None")
    mp_prec = Teuchos::null;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner method " << prec_method
		       << "." << std::endl);

  return mp_prec;
}

Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>
Stokhos::MPPreconditionerFactory::
buildPointPreconditionerFactory()
{
  std::string prec_name = 
    params->get("MP Preconditioner Type", "Ifpack");
  Teuchos::RCP<Teuchos::ParameterList> precParams = 
    Teuchos::rcp(&params->sublist("MP Preconditioner Parameters"),false);
  return Teuchos::rcp(new Stokhos::PreconditionerFactory(prec_name, 
							 precParams));
}
