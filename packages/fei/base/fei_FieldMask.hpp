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


#ifndef _fei_FieldMask_hpp_
#define _fei_FieldMask_hpp_

#include "fei_macros.hpp"

#include <vector>

namespace fei {

  /** Internal implementation class.
      A FieldMask describes the layout of the fields associated with an
      identifier.
      Each identifier Record should have a pointer to a FieldMask.
      Each FieldMask will generally be shared by a large number of identifier
      records. (i.e., the number of field-mask objects in memory will be
      small -- it will be the number of distinct fields-per-node combinations
      in the mesh. Example: If a problem has two fields (say temperature and
      displacement) and some nodes have 1 field, some have the other, and some
      have both, then there would be a total of 3 field-masks.
   */
  class FieldMask {
  public:
    /** Default Constructor */
    FieldMask();

    /** Copy Constructor */
    FieldMask(const FieldMask& fm);

    /** Construct with initial fields */
    FieldMask(int numFields,
	      const int* fieldIDs,
	      const int* fieldSizes);

    /** Destructor */
    virtual ~FieldMask();

    /** Return the mask-id for this object. */
    int getMaskID() { return( maskID_ ); }

    /** Return the mask-id that corresponds to the specified data. */
    static int calculateMaskID(int numFields,
			       const int* fieldIDs);

    /** Return the mask-id that corresponds to the specified data. */
    static int calculateMaskID(const FieldMask& fm,
			       int fieldID);

    /** Add a field-id to this object. */
    void addField(int fieldID,
		  int fieldSize);

    /** Query whether the specified field is contained in this field-mask.*/
    bool hasFieldID(int fieldID) const
    {
      for(size_t i=0; i<fieldIDs_.size(); ++i)
        if (fieldIDs_[i] == fieldID) return true;
      return false;
    }

    /** Return the number of fields in this field-mask. */
    size_t getNumFields() const { return(fieldIDs_.size()); }

    /** Return the number of global indices corresponding to the set of fields
	represented by this mask. This is sum(fieldSizes[i]).
    */
    int getNumIndices() const { return( numIndices_ ); }

    /** Set the number of global indices corresponding to the set of fields
	represented by this mask.
    */
    void setNumIndices(int numInd) { numIndices_ = numInd; }

    /** Return an array of the fields in this field-mask. */
    std::vector<int>& getFieldIDs() { return(fieldIDs_); }

    /** Return an array of the fields in this field-mask. */
    const std::vector<int>& getFieldIDs() const { return(fieldIDs_); }

    /** Return an array of the fieldSizes in this field-mask. */
    std::vector<int>& getFieldSizes() { return(fieldSizes_); }

    /** Return an array of the fieldSizes in this field-mask. */
    const std::vector<int>& getFieldSizes() const { return(fieldSizes_); }

    /** Given a field-id, return the offset of the corresponding equation-
        number in a record's list of equation-numbers.
    */
    void getFieldEqnOffset(int fieldID, int& offset) const;

    /** Test equality of field-masks. */
    bool operator==(const FieldMask& fm) const
      { return( maskID_ == fm.maskID_ ); }

    /** Test inequality of field-masks. */
    bool operator!=(const FieldMask& fm) const
      { return( maskID_ != fm.maskID_ ); }

    /** Query whether 'this' is a subset of the input argument. i.e., query
     * whether all of the fields in 'this' field-mask occur in the input
     * field-mask.
     */
    bool isSubSetOf(const FieldMask& fm) const
      {
        for(unsigned i=0; i<fieldIDs_.size(); ++i) {
          if (!fm.hasFieldID(fieldIDs_[i])) return(false);
        }
        return(true);
      }

  private:
    FieldMask& operator=(const FieldMask& src);

    int calculateMaskID();

    int maskID_;

    std::vector<int> fieldIDs_;
    std::vector<int> fieldSizes_;
    std::vector<int> fieldEqnOffsets_;

    int numFields_;
    int numIndices_;
  };

} //namespace fei

#endif // _fei_FieldMask_hpp_

