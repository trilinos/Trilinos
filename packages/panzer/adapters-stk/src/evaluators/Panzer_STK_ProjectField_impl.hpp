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

#ifndef PANZER_STK_PROJECT_FIELD_IMPL_HPP
#define PANZER_STK_PROJECT_FIELD_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_ProjectionTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

namespace panzer_stk {

template<typename EvalT,typename Traits>
ProjectField<EvalT, Traits>::
ProjectField(const std::string & inName, Teuchos::RCP<panzer::PureBasis> src,
             Teuchos::RCP<panzer::PureBasis> dst, std::string outName):
  srcBasis_(src), dstBasis_(dst)
{ 
  using panzer::Cell;
  using panzer::BASIS;

  static_assert(std::is_same<EvalT,panzer::Traits::Residual>::value);

  Teuchos::RCP<PHX::DataLayout> srcBasis_layout = srcBasis_->functional;
  Teuchos::RCP<PHX::DataLayout> dstBasis_layout = dstBasis_->functional;

  if (outName == "") outName = inName; 
  result_ = PHX::MDField<ScalarT,Cell,BASIS>(outName,dstBasis_layout);
  this->addEvaluatedField(result_);

  // This shouldn't get modified but needs to be non const 
  // because of downstream templating in intrepid2
  source_ = PHX::MDField<ScalarT,Cell,BASIS>(inName,srcBasis_layout);
  this->addNonConstDependentField(source_);

  this->setName("Project Field");

  // storage for local (to the workset) orientations
  auto maxWorksetSize = srcBasis_->functional->extent(0); // max number of cells
  local_orts_ = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orts",maxWorksetSize);

}

// **********************************************************************
template<typename EvalT,typename Traits>
void ProjectField<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData  d, 
		      PHX::FieldManager<Traits>& /* fm */)
{
  // coming from the orientations interface, this includes all orts for the process
  orientations_ = d.orientations_;
}

// **********************************************************************
template<typename EvalT,typename Traits>
void ProjectField<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  // Is there a chance workset is empty?
  if (workset.num_cells<=0) return;

  // Perform local L2 projection
#ifdef HAVE_INTREPID2_EXPERIMENTAL_NAMESPACE
  using pts = Intrepid2::Experimental::ProjectionTools<PHX::Device>;
#else
  using pts = Intrepid2::ProjectionTools<PHX::Device>;
#endif

  size_t numCells = workset.num_cells;
  const auto cell_range = std::pair<int,int>(0,numCells);
  // local_orts_ may be too large for final workset so we do the standard subview trick
  auto sub_local_orts = Kokkos::subview(local_orts_,cell_range);
  auto orts_host = Kokkos::create_mirror_view(sub_local_orts);

  // First, need to copy orientations to device
  if (orientations_ == Teuchos::null) {
    // If your bases don't require orientations, pass the default (0,0) orientation
    for (size_t i=0; i < numCells; ++i)
      orts_host(i) = Intrepid2::Orientation();
  } else {
    for (size_t i=0; i < numCells; ++i) // grab orientations for this workset
      orts_host(i) = orientations_->at(workset.cell_local_ids[i]);
  }   
  Kokkos::deep_copy(sub_local_orts,orts_host);

  // TODO BWR Revisit this... maybe we don't need pure basis upstream?
  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > dstBasis = dstBasis_->getIntrepid2Basis();
  Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > srcBasis = srcBasis_->getIntrepid2Basis();

  // Same here, need subviews
  auto sub_result = Kokkos::subview(result_.get_view(),cell_range,Kokkos::ALL());
  auto sub_source = Kokkos::subview(source_.get_view(),cell_range,Kokkos::ALL());

  pts::projectField(sub_result,dstBasis.get(),
                    sub_source,srcBasis.get(),sub_local_orts);

}

}

#endif