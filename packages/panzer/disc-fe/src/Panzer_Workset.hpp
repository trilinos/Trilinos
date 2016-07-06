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


#ifndef PANZER_WORKSET_HPP
#define PANZER_WORKSET_HPP

#include <cstddef>
#include <vector>
#include <map>
#include <iostream>

#include "Panzer_Dimension.hpp"
#include "Panzer_BasisValues2.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_Dimension.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

  class LinearObjContainer;

  /** This is used within the workset to make edge based (DG like) assembly
    * an easier task. This basically allows seperation of the workset abstraction
    * from how it is accessed.
    */
  struct WorksetDetails {
    typedef PHX::MDField<double,Cell,NODE,Dim> CellCoordArray;

    Kokkos::View<const int*,PHX::Device> cell_local_ids_k;
    std::vector<std::size_t> cell_local_ids;
    CellCoordArray cell_vertex_coordinates;
    std::string block_id;

    int subcell_index; //! If workset corresponds to a sub cell, what is the index?

    //! Value correspondes to integration order.  Use the offest for indexing.
    Teuchos::RCP< std::vector<int> > ir_degrees;
    
    std::vector<Teuchos::RCP<panzer::IntegrationValues2<double> > > int_rules;
    
    //! Value corresponds to basis type.  Use the offest for indexing.
    Teuchos::RCP< std::vector<std::string> > basis_names;

    //! Static basis function data, key is basis name, value is index in the static_bases vector
    std::vector<Teuchos::RCP< panzer::BasisValues2<double> > > bases;
  };

  /** This is the main workset object. Not that it inherits from WorksetDetails, this
    * is to maintain backwards compatibility in the use of the workset object. The addition
    * of a details vector supports things like DG based assembly.
    */
  struct Workset : public WorksetDetails {
    Workset() {}

    std::size_t num_cells;
    int subcell_dim; //! If workset corresponds to a sub cell, what is the dimension?
    
    double alpha;
    double beta;
    double time;
    std::vector<double> gather_seeds; // generic gather seeds
    bool evaluate_transient_terms;

    //! other contains details about the side-sharing elements on the other side
    //! of the interface. If Teuchos::nonnull(other), then Workset contains two
    //! WorksetDetails: itself, and other.
    Teuchos::RCP<WorksetDetails> other;

    //! op(0) return *this; op(1) returns *other.
    WorksetDetails& operator()(const int i) {
      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
      return i == 0 ? static_cast<WorksetDetails&>(*this) : *other;
    }
    //! const accessor.
    const WorksetDetails& operator()(const int i) const {
      TEUCHOS_ASSERT(i == 0 || (i == 1 && Teuchos::nonnull(other)));
      return i == 0 ? static_cast<const WorksetDetails&>(*this) : *other;
    }
    //! Convenience wrapper to operator() for pointer access.
    WorksetDetails& details(const int i) { return operator()(i); }
    const WorksetDetails& details(const int i) const { return operator()(i); }
    //! Return the number of WorksetDetails this Workset holds.
    size_t numDetails() const { return Teuchos::nonnull(other) ? 2 : 1; }
  };

  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w);

  /** This accessor may be used by an evaluator to abstract Details
    * Index. "Details Index" may be in the evaluator's constructor
    * ParameterList. If it is, then this accessor's constructor reads it and
    * checks for error. Regardless of whether it is, operator() then wraps
    * workset's access to WorksetDetails fields to provide the correct ones.
    */
  class WorksetDetailsAccessor {
  public:
    //! Default value is 0, which is backwards compatible.
    WorksetDetailsAccessor() : details_index_(0) {}
    //! An evaluator builder sets the details index.
    void setDetailsIndex(const int di) { details_index_ = di; }
    //! Get the details index. Generally, only a gather evaluator needs to know
    //! its details index.
    int getDetailsIndex() const { return details_index_; }
    //! Workset wrapper to extract the correct details. Example: wda(workset).bases[i].
    WorksetDetails& operator()(Workset& workset) const {
      return workset(details_index_);
    }
    //! const accessor.
    const WorksetDetails& operator()(const Workset& workset) const {
      return workset(details_index_);
    }
  private:
    int details_index_;
  };

} // namespace panzer

#endif
