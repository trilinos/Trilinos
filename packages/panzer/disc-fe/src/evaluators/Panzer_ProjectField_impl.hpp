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

#ifndef PANZER_PROJECT_FIELD_IMPL_HPP
#define PANZER_PROJECT_FIELD_IMPL_HPP

// TODO BWR check this file and decl for unnecessary includes

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::ProjectField<EvalT, Traits>::
ProjectField(const Teuchos::ParameterList& p)
{ 
  // TODO HERE could initialize differently, without param list!!
  field_name_ = (p.get< std::string >("Field Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Source Basis"))
    srcBasis_ = p.get< Teuchos::RCP<PureBasis> >("Source Basis");
  else
    srcBasis_ = p.get< Teuchos::RCP<const PureBasis> >("Source Basis");

  if(p.isType< Teuchos::RCP<PureBasis> >("Destination Basis"))
    dstBasis_ = p.get< Teuchos::RCP<PureBasis> >("Destination Basis");
  else
    dstBasis_ = p.get< Teuchos::RCP<const PureBasis> >("Destination Basis");

  Teuchos::RCP<PHX::DataLayout> srcBasis_layout = srcBasis_->functional;
  Teuchos::RCP<PHX::DataLayout> dstBasis_layout = dstBasis_->functional;

  // TODO OK to have same name?
  result_ = PHX::MDField<double,Cell,BASIS>(field_name_,dstBasis_layout);
  this->addEvaluatedField(result_);

  // This shouldn't get modified but needs to be non const 
  // because of downstream templating in intrepid2
  source_ = PHX::MDField<double,Cell,BASIS>(field_name_,srcBasis_layout);
  this->addNonConstDependentField(source_);

  this->setName("Project Field");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectField<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  orientations_ = d.orientations_;
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::ProjectField<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  // TODO is there a chance workset is empty?
  if (workset.num_cells<=0) return;

  // Perform local L2 projection

  typedef Intrepid2::Experimental::ProjectionTools<PHX::Device> pts;

  // First, need to copy orientations to device
  auto orts = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts",orientations_->size());
  auto orts_host = Kokkos::create_mirror_view(orts);
  for (size_t i=0; i < orts_host.extent(0); ++i)
    orts_host(i) = orientations_->at(i);
  Kokkos::deep_copy(orts,orts_host);

  // TODO BWR maybe we dont need pure basis upstream?
  // TODO BWR ok to make this all double here and for MD fields above?
  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > dstBasis = dstBasis_->getIntrepid2Basis();
  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > srcBasis = srcBasis_->getIntrepid2Basis();

  pts::projectField(result_.get_view(),dstBasis.get(),
                    source_.get_view(),srcBasis.get(),orts);

}

#endif
