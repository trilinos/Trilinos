/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

