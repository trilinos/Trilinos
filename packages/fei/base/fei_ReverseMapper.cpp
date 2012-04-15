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


#include <fei_ReverseMapper.hpp>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <fei_VectorSpace.hpp>
#include <snl_fei_RecordCollection.hpp>

namespace fei {

ReverseMapper::ReverseMapper(const VectorSpace& vspace)
 : eqnmap_()
{
  std::vector<int> idTypes;
  vspace.getIDTypes(idTypes);

  const std::vector<int>& eqnNumbers = vspace.getEqnNumbers();

  for(size_t idt=0; idt<idTypes.size(); ++idt) {
    const snl_fei::RecordCollection* recordcollection = NULL;
    int err = vspace.getRecordCollection(idTypes[idt], recordcollection);
    if (err != 0) {
      throw std::runtime_error("fei::ReverseMapper ERROR, failed to retrieve record-collection.");
    }

    const std::vector<fei::Record<int> >&
        records = recordcollection->getRecords();

    for(size_t i=0; i<records.size(); ++i) {
      const fei::Record<int>* record = &records[i];

      const fei::FieldMask* fm = record->getFieldMask();
      const std::vector<int>& fieldIDs = fm->getFieldIDs();
      const std::vector<int>& fieldSizes = fm->getFieldSizes();

      int offsetIntoEqnNumbers = record->getOffsetIntoEqnNumbers();

      for(size_t i=0; i<fieldIDs.size(); ++i) {
        int offset2 = 0;
        fm->getFieldEqnOffset(fieldIDs[i], offset2);

        EqnRecord erec;
        erec.IDType = idTypes[idt];
        erec.ID = record->getID();
        erec.fieldID = fieldIDs[i];

        for(int j=0; j<fieldSizes[i]; ++j) {
          erec.offset = j;
          erec.global_eqn = eqnNumbers[offsetIntoEqnNumbers+offset2+j];
          eqnmap_.insert(std::make_pair(erec.global_eqn, erec));
        }
      }
    }
  }
}

ReverseMapper::~ReverseMapper()
{
}

EqnRecord ReverseMapper::getEqnRecord(int global_eqn, int option) const
{
  std::map<int,EqnRecord>::const_iterator
    iter = eqnmap_.find(global_eqn);

  if (iter == eqnmap_.end()) {
    if (option == 1) {
      EqnRecord erec;
      erec.IDType = -1;
      erec.ID = -1;
      erec.fieldID = -1;
      erec.offset = -1;
      erec.global_eqn = -1;
      return erec;
    }

    std::ostringstream osstr;
    osstr << "fei::ReverseMapper::getEqnRecord ERROR, global_eqn="<<global_eqn
       << " not found.";
    throw std::runtime_error(osstr.str());
  }

  return iter->second;
}

}//namespace fei

