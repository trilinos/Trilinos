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
#include "Teuchos_ArrayRCP.hpp"
#include "Panzer_Dimension.hpp"
#include "Shards_Array.hpp"
#include "Panzer_BasisValues.hpp"
#include "Panzer_IntegrationValues.hpp"

#include "Intrepid_FieldContainer.hpp"

namespace panzer {

  class LinearObjContainer;

  /** This is used within the workset to make edge based (DG like) assembly
    * an easier task. This basically allows seperation of the workset abstraction
    * from how it is accessed.
    */
  struct WorksetDetails {
    std::vector<std::size_t> cell_local_ids;
    Intrepid::FieldContainer<double> cell_vertex_coordinates;
    std::string block_id;

    int subcell_index; //! If workset corresponds to a sub cell, what is the index?

    //! Value correspondes to integration order.  Use the offest for indexing.
    Teuchos::RCP< std::vector<int> > ir_degrees;
    
    std::vector<Teuchos::RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > > int_rules;
    
    //! Value corresponds to basis type.  Use the offest for indexing.
    Teuchos::RCP< std::vector<std::string> > basis_names;

    //! Static basis function data, key is basis name, value is index in the static_bases vector
    std::vector<Teuchos::RCP< panzer::BasisValues<double,Intrepid::FieldContainer<double> > > > bases;
  };

  /** This is the main workset object. Not that it inherits from WorksetDetails, this
    * is to maintain backwards compatibility in the use of the workset object. The addition
    * of a details vector supports things like DG based assembly.
    */
  struct Workset : public WorksetDetails {
    
    std::size_t num_cells;
    int subcell_dim; //! If workset corresponds to a sub cell, what is the dimension?

    std::vector<Teuchos::RCP<WorksetDetails> > details; // note that (*this) is set to index [0]
                                                        // when using panzers workset construction
    
    double alpha;
    double beta;
    double time;
    std::vector<double> gather_seeds; // generic gather seeds
    bool evaluate_transient_terms;
  };

  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w);

} // namespace panzer

#endif
