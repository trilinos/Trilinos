/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_Dof_hpp_
#define _fei_Dof_hpp_

#include <fei_macros.hpp>

namespace fei {

/** Dof - mesh-degree-of-freedom.
 *
 * A mesh-dof is the triple (rank, id, field).
 * - Rank is an integer type used to label different 'kinds' of mesh-dofs,
 *   e.g. node vs face, etc.
 * - Id identifies a particular instance of a rank, e.g., node 97.
 * - Field is an integer type used to label solution fields such as
 *   temperature or displacement, etc.
 *
 * Thus if the user chooses to give nodes a rank of 0 and the temperature
 * field a label of 4, then the temperature at node 97 can be represented
 * as the dof (0, 97, 4).
 *
 * Notes:
 *
 * 1. The Dof class is templated on the two integer types LocalOrdinal
 *    and GlobalOrdinal. Ranks and Fields have type LocalOrdinal, while
 *    ids have type GlobalOrdinal. The distinction is somewhat arbitrary,
 *    but the assumption is that there may be billions (or more?) ids, but
 *    probably far fewer distinct ranks or fields. So in extreme cases a
 *    user may wish to use a different type (larger) for ids than for the
 *    rank or field.
 */
template<class LocalOrdinal, class GlobalOrdinal>
class Dof {
 public:
  /** constructor */
  Dof(LocalOrdinal rank, GlobalOrdinal id, LocalOrdinal field)
   : m_rank(rank), m_id(id), m_field(field) {}

  /** destructor */
  ~Dof(){}

  LocalOrdinal rank() const { return m_rank; }
  GlobalOrdinal id() const { return m_id; }
  LocalOrdinal field() const { return m_field; }

 private:
  LocalOrdinal m_rank;
  GlobalOrdinal m_id;
  LocalOrdinal m_field;
};//class Dof

/** Less operator which will order Dofs by rank first, then id then field.
 */
template<class LocalOrdinal, class GlobalOrdinal>
struct less_rank_id_field {
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>& dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>& dof2) const
  {
    if (dof1.rank()==dof2.rank()) {
      if (dof1.id() == dof2.id()) return dof1.field() < dof2.field();
      else return dof1.id() < dof2.id();
    }
    else {
      return dof1.rank() < dof2.rank();
    }
  }

  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>* dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>* dof2) const
  {
    if (dof1->rank()==dof2->rank()) {
      if (dof1->id() == dof2->id()) return dof1->field() < dof2->field();
      else return dof1->id() < dof2->id();
    }
    else {
      return dof1->rank() < dof2->rank();
    }
  }
};

/** Less operator which will order Dofs by field first, then rank then id.
 */
template<class LocalOrdinal, class GlobalOrdinal>
struct less_field_rank_id {
  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>& dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>& dof2) const
  {
    if (dof1.field()==dof2.field()) {
      if (dof1.rank() == dof2.rank()) return dof1.id() < dof2.id();
      else return dof1.rank() < dof2.rank();
    }
    else {
      return dof1.field() < dof2.field();
    }
  }

  bool operator()(const Dof<LocalOrdinal,GlobalOrdinal>* dof1,
                 const Dof<LocalOrdinal,GlobalOrdinal>* dof2) const
  {
    if (dof1->field()==dof2->field()) {
      if (dof1->rank() == dof2->rank()) return dof1->id() < dof2->id();
      else return dof1->rank() < dof2->rank();
    }
    else {
      return dof1->field() < dof2->field();
    }
  }
};

}//namespace fei
#endif

