/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

