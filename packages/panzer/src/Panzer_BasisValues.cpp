// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_config.hpp"
#include "Panzer_BasisValues.hpp"
#include "Panzer_BasisValues_impl.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_Traits.hpp"


namespace panzer {

template<typename Scalar,typename Array>
template <typename ArrayFactory>
void panzer::BasisValues<Scalar,Array>::
setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& layout,
            const ArrayFactory & af)
{
  basis_layout = layout;
  Teuchos::RCP<const panzer::PureBasis> basisDesc = layout->getBasis();

  // for convience pull out basis and quadrature information
  int num_quad = layout->getNumPoints();
  int dim      = basisDesc->getDimension();
  int card     = basisDesc->getCardinality();
  int numcells = basisDesc->getNumCells();
  panzer::PureBasis::EElementSpace elmtspace = basisDesc->getElementSpace();
  Teuchos::RCP<const shards::CellTopology> cellTopo = basisDesc->getCellTopology();
  
  intrepid_basis = basisDesc->getIntrepidBasis<Scalar,Array>();
  
  // allocate field containers
  // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec

  // compute basis fields
  if(elmtspace==panzer::PureBasis::HGRAD) {
     // HGRAD is a nodal field

     // build values
     ///////////////////////////////////////////////

     basis_ref = af.template buildArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
     basis = af.template buildArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

     weighted_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

     // build gradients
     ///////////////////////////////////////////////

     grad_basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("grad_basis_ref",card,num_quad,dim); // F, P, D
     grad_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("grad_basis",numcells,card,num_quad,dim);

     weighted_grad_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_grad_basis",numcells,card,num_quad,dim);

     // build curl
     ///////////////////////////////////////////////

     // nothing - HGRAD does not support CURL operation
  }
  else if(elmtspace==panzer::PureBasis::HCURL) {
     // HCURL is a vector field

     // build values
     ///////////////////////////////////////////////

     basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
     basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

     weighted_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - HCURL does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     if(dim==2) {
        // curl of HCURL basis is not dimension dependent
        curl_basis_ref = af.template buildArray<Scalar,BASIS,IP>("curl_basis_ref",card,num_quad); // F, P
        curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("curl_basis",numcells,card,num_quad);

        weighted_curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP>("weighted_curl_basis",numcells,card,num_quad);
     }
     else if(dim==3){
        curl_basis_ref = af.template buildArray<Scalar,BASIS,IP,Dim>("curl_basis_ref",card,num_quad,dim); // F, P, D
        curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("curl_basis",numcells,card,num_quad,dim);

        weighted_curl_basis = af.template buildArray<Scalar,Cell,BASIS,IP,Dim>("weighted_curl_basis",numcells,card,num_quad,dim);
     }
     else { TEUCHOS_ASSERT(false); } // what do I do with 1D?
  }
  else { TEUCHOS_ASSERT(false); }

  basis_coordinates_ref = af.template buildArray<Scalar,BASIS,Dim>("basis_coordinates_ref",card,dim);
  basis_coordinates = af.template buildArray<Scalar,Cell,BASIS,Dim>("basis_coordinates",numcells,card,dim);
}

// do some explicit instantiation so things run faster.

#define BASIS_VALUES_INSTANTIATION(SCALAR) \
template class BasisValues<SCALAR,Intrepid::FieldContainer<SCALAR> >;\
template class BasisValues<SCALAR,PHX::MDField<SCALAR> >;\
template \
void BasisValues<SCALAR,Intrepid::FieldContainer<SCALAR> >::setupArrays<IntrepidFieldContainerFactory>( \
                   const Teuchos::RCP<const panzer::BasisIRLayout>& basis, \
                   const IntrepidFieldContainerFactory & af); \
template \
void BasisValues<SCALAR,PHX::MDField<SCALAR> >::setupArrays<MDFieldArrayFactory>( \
                   const Teuchos::RCP<const panzer::BasisIRLayout>& basis, \
                   const MDFieldArrayFactory & af);

BASIS_VALUES_INSTANTIATION(panzer::Traits::RealType)
BASIS_VALUES_INSTANTIATION(panzer::Traits::FadType)
#ifdef HAVE_STOKHOS
  BASIS_VALUES_INSTANTIATION(panzer::Traits::SGType)
  BASIS_VALUES_INSTANTIATION(panzer::Traits::SGFadType)
#endif

} // namespace panzer
