#ifndef _fei_BCRecord_hpp_
#define _fei_BCRecord_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"

/**  BCRecord is a boundary condition specification for one field on one node.
     This class defines comparison operators (<,>,==,!=) which only compare the
     nodeIDs and fieldIDs, ignoring alpha/beta/gamma. In other words, these
     comparison operators are to be used for creating lists of BCRecords that
     specify BCs for distinct nodeID/fieldID pairs. Obviously equality occurs
     when both nodeID and fieldID are the same. The less-than and greater-than
     operators assume that BCRecords are sorted by nodeID first, and fieldID
     second. In other words, a BCRecord on nodeID 1, fieldID 2 will be deemed
     "less than" a BCRecord on nodeID 2, fieldID 1.
 */
class BCRecord {
 public:
  /** Default constructor. */
  BCRecord();

  virtual ~BCRecord() {}

  /** Set the incoming data onto this BCRecord. NOTE: the coefPtr is a pointer
      to an array to be used for this BCRecord's coefs. This BCRecord will keep
      this pointer, and refer to the memory pointed to by it. So don't destroy
      this memory before destroying this BCRecord object. The memory pointed to
      by coefPtr will be assumed to be of length 3*fieldSize.
  */
  int init(GlobalID nodeID, int fieldID, int fieldSize, double* coefPtr);

  int getIDType() const {return(idType_);}
  void setIDType(int idtype) {idType_ = idtype;}

  GlobalID getNodeID() const {return(nodeID_);}
  int getFieldID() const {return(fieldID_);}
  int getFieldSize() const {return(fieldSize_);}

  void setAlpha(const double* alpha) { setDoubleList(BC_ALPHA, alpha); }
  void setBeta(const double* beta)   { setDoubleList(BC_BETA,  beta ); }
  void setGamma(const double* gamma) { setDoubleList(BC_GAMMA, gamma); }

  const double *pointerToAlpha() const {return( &(coef_[BC_ALPHA * fieldSize_]) );}
  const double *pointerToBeta()  const {return( &(coef_[BC_BETA  * fieldSize_]) );}
  const double *pointerToGamma() const {return( &(coef_[BC_GAMMA * fieldSize_]) );}

  bool operator==(const BCRecord& bcrec) const
    {
      if ( (nodeID_ == bcrec.nodeID_) && (fieldID_ == bcrec.fieldID_) )
	return(true);
      return(false);
    }

  bool operator!=(const BCRecord& bcrec) const
    {
      if ( (nodeID_ == bcrec.nodeID_) && (fieldID_ == bcrec.fieldID_) )
	return(false);
      return(true);
    }

  bool operator<(const BCRecord& bcrec) const
    {
      if (nodeID_ < bcrec.nodeID_) return(true);
      if (nodeID_ > bcrec.nodeID_) return(false);

      //if we get to here nodeIDs must be equal, so now check fieldIDs.

      if (fieldID_ < bcrec.fieldID_) return(true);

      //now we know that either the fieldIDs are equal, or bcrec's fieldID is
      //less than ours. In either case, "<" is false.
  
      return(false);
    }

  bool operator>(const BCRecord& bcrec) const
    {
      if (nodeID_ > bcrec.nodeID_) return(true);
      if (nodeID_ < bcrec.nodeID_) return(false);

      //if we get to here nodeIDs must be equal, so now check fieldIDs.

      if (fieldID_ > bcrec.fieldID_) return(true);

      //now we know that either the fieldIDs are equal, or bcrec's fieldID is
      //greater than ours. In either case, ">" is false.
  
      return(false);
    }

 private:
  BCRecord(const BCRecord& src);

  BCRecord& operator=(const BCRecord& src);

  void setDoubleList(int bcComponent, const double* input);

   int idType_;
   GlobalID nodeID_;
   int fieldID_;                       // cached field ID
   int fieldSize_;                     // cached field cardinality

   enum { BC_ALPHA = 0, BC_BETA = 1, BC_GAMMA = 2 };

   double *coef_;     // alpha/beta/gamma coefficients
};
 
#endif
