/*
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
*/

#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "Belos_TpetraAdapter_UQ_PCE.hpp"
#include "Stokhos_MueLu_UQ_PCE.hpp"

#include <Kokkos_Core.hpp>
#include <HexElement.hpp>
#include <fenl_functors_pce.hpp>
#include <fenl_impl.hpp>

namespace Kokkos {
namespace Example {
namespace FENL {

typedef Stokhos::DynamicStorage<int,double,Cuda> Storage;
typedef Sacado::UQ::PCE<Storage> Scalar;
typedef ElementComputationKLCoefficient<Scalar,double,Cuda> KL_Vector;
typedef ElementComputationKLCoefficient<double,double,Cuda> KL_Scalar;

INST_FENL( Scalar , Cuda , BoxElemPart::ElemLinear ,
           KL_Vector , TrivialManufacturedSolution )
INST_FENL( Scalar , Cuda , BoxElemPart::ElemQuadratic ,
           KL_Vector , TrivialManufacturedSolution )
INST_KL( Scalar , double , Cuda )

INST_FENL( double , Cuda , BoxElemPart::ElemLinear ,
           KL_Scalar , TrivialManufacturedSolution )
INST_FENL( double , Cuda , BoxElemPart::ElemQuadratic ,
           KL_Scalar , TrivialManufacturedSolution )
INST_KL( double , double , Cuda )

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos */
