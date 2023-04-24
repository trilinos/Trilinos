// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef SOLO_SIDE_ID_GENERATOR_HPP_
#define SOLO_SIDE_ID_GENERATOR_HPP_

#include <stk_util/util/ReportHandler.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {
namespace impl {

class SoloSideIdGenerator
{
public:
    SoloSideIdGenerator(int nProcs, int proc, uint64_t maxSideId)
      : numProcs(nProcs), pseudoOrdinal(6)
    {
        calculate_max_pseudo_element(maxSideId);
        pseudoElement = proc + 1;
    }

    void use_fmwk_id_type_to_determine_largest_valid_id()
    {
        calculate_max_pseudo_element(std::numeric_limits<int>::max());
    }

    stk::mesh::EntityId get_solo_side_id()
    {
        stk::mesh::EntityId id = get_solo_side_id_using_formula(pseudoElement, pseudoOrdinal);
        STK_ThrowRequireMsg(pseudoElement <= max_pseudo_element(), "Exhausted solo side ids. Please report to sierra-help@sandia.gov");
        incrementPseudoElementAndOrdinal();
        return id;
    }

    uint64_t max_pseudo_element() const  { return m_maxPseudoElement; }

protected:
    stk::mesh::EntityId get_solo_side_id_using_formula(unsigned elementId, unsigned sideOrdinal)
    {
        STK_ThrowRequireMsg(elementId <= m_maxPseudoElement, "Exhausted solo side ids for this processor. elementId= " << elementId << ", m_maxPseudoElement= " 
            << m_maxPseudoElement << ". Please report to sierra-help@sandia.gov");

        //this is the side-id formula used by IO. the "+1" is because IO always uses one-based side ordinals
        return 10*elementId + sideOrdinal + 1;
    }

private:
    void calculate_max_pseudo_element(uint64_t maxSideId)
    {
        STK_ThrowRequire(maxSideId > 10);
        m_maxPseudoElement = (maxSideId - 10) / 10;
    }
    void incrementPseudoElementAndOrdinal()
    {
        pseudoOrdinal++;
        if( pseudoOrdinal % 10 == 0)
        {
            pseudoOrdinal = 6;
            pseudoElement += numProcs;
        }
    }
    int numProcs;
    uint64_t m_maxPseudoElement;
    unsigned pseudoOrdinal;
    stk::mesh::EntityId pseudoElement;
};

}
}
}

#endif /* SOLO_SIDE_ID_GENERATOR_HPP_ */
