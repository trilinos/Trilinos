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

#include "Panzer_Workset.hpp"

namespace panzer {
 
  std::ostream& operator<<(std::ostream& os, const panzer::Workset& w)
  {
    using std::endl;

    os << "Workset" << endl;
    os << "  block_id=" << w.block_id << endl;
    os << "  num_cells:" << w.num_cells << endl;
    os << "  cell_local_ids (size=" << w.cell_local_ids.size() << ")" << endl;
    os << "  subcell_dim = " << w.subcell_dim << endl;
    os << "  subcell_index = " << w.subcell_index << endl;

    os << "  ir_degrees: " << endl;
    for (std::vector<int>::const_iterator ir = w.ir_degrees->begin();
	 ir != w.ir_degrees->end(); ++ir)
      os << "    " << *ir << std::endl;

    std::vector<int>::const_iterator ir = w.ir_degrees->begin();
    for (std::vector<Teuchos::RCP<panzer::IntegrationValues<double,Intrepid::FieldContainer<double> > > >::const_iterator irv = w.int_rules.begin();  irv != w.int_rules.end(); ++irv,++ir) {

      os << "  IR Values (Degree=" << *ir << "):" << endl;

      os << "    cub_points:" << endl;
      os << (*irv)->cub_points << endl;

      os << "    side_cub_points:" << endl;
      os << (*irv)->side_cub_points << endl;

      os << "    cub_weights:" << endl;
      os << (*irv)->cub_weights << endl;

      os << "    node_coordinates:" << endl;
      os << (*irv)->node_coordinates << endl;

      os << "    jac:" << endl;
      os << (*irv)->jac << endl;

      os << "    jac_inv:" << endl;
      os << (*irv)->jac_inv << endl;

      os << "    jac_det:" << endl;
      os << (*irv)->jac_det << endl;

      os << "    weighted_measure:" << endl;
      os << (*irv)->weighted_measure << endl;

      os << "    covarient:" << endl;
      os << (*irv)->covarient << endl;

      os << "    contravarient:" << endl;
      os << (*irv)->contravarient << endl;

      os << "    norm_contravarient:" << endl;
      os << (*irv)->norm_contravarient << endl;

      os << "    ip_coordinates:" << endl;
      os << (*irv)->ip_coordinates << endl;

      os << "    int_rule->getName():" << (*irv)->int_rule->getName() << endl;
    }


    os << "  basis_names: " << endl;
    for (std::vector<std::string>::const_iterator b = w.basis_names->begin();
	 b != w.basis_names->end(); ++b)
      os << "    " << *b << std::endl;

    std::vector<std::string>::const_iterator b = w.basis_names->begin();
    for (std::vector<Teuchos::RCP< panzer::BasisValues<double,Intrepid::FieldContainer<double> > > >::const_iterator bv = w.bases.begin(); bv != w.bases.end(); ++bv,++b) {

      os << "  Basis Values (basis_name=" << *b << "):" << endl;
      
      os << "    basis_ref:" << endl;
      os << (*bv)->basis_ref << endl;

      os << "    basis:" << endl;
      os << (*bv)->basis << endl;

      os << "    grad_basis_ref:" << endl;
      os << (*bv)->grad_basis_ref << endl;

      os << "    grad_basis:" << endl;
      os << (*bv)->grad_basis << endl;

      os << "    curl_basis_ref:" << endl;
      os << (*bv)->curl_basis_ref << endl;

      os << "    curl_basis:" << endl;
      os << (*bv)->curl_basis << endl;

      os << "    basis_coordinates_ref:" << endl;
      os << (*bv)->basis_coordinates_ref << endl;

      os << "    basis_coordinates:" << endl;
      os << (*bv)->basis_coordinates << endl;

      os << "    basis_layout->name():" << (*bv)->basis_layout->name() << endl;
    }

  

    return os;
  }

}
