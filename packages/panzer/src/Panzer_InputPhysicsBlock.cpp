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


#include <vector>
#include <string>
#include <iostream>
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_InputEquationSet.hpp"

std::ostream& 
panzer::operator<<(std::ostream& os, const panzer::InputPhysicsBlock& i) {
  
  os << "Physics Block: " << i.physics_block_id 
     << ", Num Equation Sets: " << i.eq_sets.size() << std::endl;
  for (std::size_t k=0; k < i.eq_sets.size(); ++k) {
    os << "  Set " << k << std::endl;
    os << "    name = " << i.eq_sets[k].name << std::endl;
    os << "    basis = " << i.eq_sets[k].basis << std::endl;
    os << "    int order = " << i.eq_sets[k].integration_order << std::endl;
    os << "    model = " << i.eq_sets[k].model_id << std::endl;
    os << "    prefix = " << i.eq_sets[k].prefix << std::endl;
    os << "    parameterlist: " << std::endl;

    Teuchos::ParameterList::PrintOptions po;
    po.showTypes(true);
    po.indent(6);
    i.eq_sets[k].params.print(os,po);
  }
  return os;
}
