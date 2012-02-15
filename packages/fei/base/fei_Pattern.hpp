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


#ifndef _fei_Pattern_hpp_
#define _fei_Pattern_hpp_

#include "fei_macros.hpp"

#include <vector>

namespace snl_fei {
  class RecordCollection;
}

namespace fei {

  /** Stencil-like pattern definition/description.
      Describes the layout of a set of field-identifiers associated with a
      set of identifiers and identifier-types.

      Example: Can be used to describe the layout of nodes with associated
      fields on an element ('element' as in finite-elements).
   */
  class Pattern {
  public:
    /** enumeration of different pattern-types */
    enum PatternType { NO_FIELD, SIMPLE, SINGLE_IDTYPE, GENERAL };

    /** Constructor, Pattern::PatternType == NO_FIELD */
    Pattern(int numIDs, int idType, snl_fei::RecordCollection* records);

    /** Constructor, Pattern::PatternType == SIMPLE
	There is only one id-type, and only one field.
     */
    Pattern(int numIDs, int idType, snl_fei::RecordCollection* records,
	    int fieldID, int fieldSize);

    /** Constructor, Pattern::PatternType == SINGLE_IDTYPE
	There is only one id-type, but there may be multiple fields per id.
     */
    Pattern(int numIDs, int idType, snl_fei::RecordCollection* records,
	    const int* numFieldsPerID,
	    const int* fieldIDs, const int* fieldSizes);

    /** Constructor, Pattern::PatternType == GENERAL
      There may be multiple id-types as well as multiple fields-per-id.
     */
    Pattern(int numIDs, const int* idTypes, snl_fei::RecordCollection*const* records,
	    const int* numFieldsPerID,
	    const int* fieldIDs, const int* fieldSizes);

    virtual ~Pattern();

    /** Return pattern-type-identifying enum
     */
    PatternType getPatternType() const { return( type_ ); }

    /** Return the number of identifiers described by this pattern. */
    int getNumIDs() const { return( numIDs_ ); }

    /** Return pointer to list of length getNumIDs() */
    const int* getIDTypes() const { return( idTypes_ ); }

    /** Return pointer to list of length getNumIDs() */
    snl_fei::RecordCollection*const* getRecordCollections() const { return &recordCollections_[0]; }

    /** Return list of length getNumIDs() */
    const int* getNumFieldsPerID() const { return( numFieldsPerID_ ); }

    /** Return list of length getTotalNumFields() */
    const int* getFieldIDs() const { return( fieldIDs_ ); }

    /** Return list of length getNumIDs() */
    const int* getNumIndicesPerID() const
    {
      return( numIndicesPerID_ );
    }

    /** total-num-fields = sum(numFieldsPerID) */
    int getTotalNumFields() const { return( totalNumFields_ ); }

    /** Return the total number of scalar indices represented by this pattern.
     This is the number of identifiers if no fields are associated with them,
     otherwise it is the sum of the field-sizes of the fields associated with
     the identifiers.
    */
    int getNumIndices() const { return( numIndices_ ); }

    /** return true if the 'rhs' pattern is the same as 'this' pattern.
     */
    bool operator==(const Pattern& rhs) const;

    /** return true if the 'rhs' pattern is different than 'this' pattern.
     */
    bool operator!=(const Pattern& rhs) const;

  private:
    PatternType type_;
    int numIDs_;
    int totalNumFields_;
    int numIndices_;
    std::vector<int> data_;
    std::vector<snl_fei::RecordCollection*> recordCollections_;

    const int* idTypes_;
    const int* numFieldsPerID_;
    const int* fieldIDs_;
    const int* numIndicesPerID_;
  };

} //namespace fei

#endif // _fei_Pattern_hpp_

