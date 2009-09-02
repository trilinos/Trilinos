/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_PointBlockMap_hpp_
#define _snl_fei_PointBlockMap_hpp_

#include <fei_macros.hpp>

#include <map>

namespace snl_fei {
  /** Stores mappings between point-entry equations and block-entry equations.
      This class, once initialized, can answer queries such as:<br>
      Given a point-equation, return the corresponding block-equation.<br>
      or, Given a block-equation, return the corresponding size (number of
      point-equations contained in the block).
  */
  class PointBlockMap {
  public:
    /** Constructor */
    PointBlockMap();

    /** Destructor */
    virtual ~PointBlockMap();

   /** Indicate to this object that point-equations are equal to
       block-equations which means that all block-equation sizes are 1 and
       all mappings are trivial. This will cause the internal mapping objects
       to be emptied, to save memory, and all lookup queries will be instant.
   */
   void setPtEqualBlk();

   /** query whether block-size == 1, i.e., "point-equals-block"
   */
   bool ptEqualBlk() { return( ptEqualBlk_ ); }

   /** Establish the mapping blkEqn => ptEqn. I.e., insert ptEqn and blkEqn 
       into internal structures, and establish a correspondence between them.
       This function returns an error if 'setPtEqualBlk' has been called and
       ptEqn != blkEqn.
   */
   int setEqn(int ptEqn, int blkEqn);

   /** Establish the mapping blkEqn => ptEqn. I.e., insert ptEqn and blkEqn 
       into internal structures, and establish a correspondence between them,
       also setting the associated block-size.
       This function returns an error if 'setPtEqualBlk' has been called and
       ptEqn != blkEqn.
   */
   int setEqn(int ptEqn, int blkEqn, int blkSize);

   /** Store the specified size corresponding to the specified block-equation.
    Note that blkEqn must already have appeared in a call to 'setEqn', in order
   to establish it in a mapping. Otherwise, an error (-1) is returned.
   */
   int setBlkEqnSize(int blkEqn, int size);

   /** Given a blkEqn, return the corresponding size. */
   int getBlkEqnSize(int blkEqn);

   /** Query the maximum block-equation size
    */
   int getMaxBlkEqnSize() { return( maxSize_ ); }

   /** Set the maximum block-equation size
    */
   void setMaxBlkEqnSize(int sz) { maxSize_ = sz; }

   /** Given a point-equation, return the corresponding block-equation.
       If eqn does not correspond to a block-equation, then -1 is
       returned.
       A return-value of -1 always indicates a not-found or not-valid error.
   */
   int eqnToBlkEqn(int eqn) const;

   /** Given a block-equation, return the corresponding point-equation (the
       first of the point-equations that correspond to that block-equation).
       If blkEqn does not correspond to a point-equation, then -1 is
       returned.
   */
   int blkEqnToPtEqn(int blkEqn) const;

   /** Given a block-equation, query for the corresponding point-equation (the
       first of the corresponding point-equations) and the block-size, or
       number of corresponding point-equations.
       @return -1 if block-equation not found, 0 if no error
   */
   int getBlkEqnInfo(int blkEqn, int& ptEqn, int& blkSize);

   /** Given a point-equation, query for the corresponding block-equation and
       the offset of this point-equation into the block-equation.
       @return -1 if point-equation not found, 0 if no error
   */
   int getPtEqnInfo(int ptEqn, int& blkEqn, int& blkOffset);

   /** Given a blkEqn/ptEqn pair, return the offset of the ptEqn into the
       blkEqn. Example: if blkEqn 23 corresponds to ptEqn 69, and the arguments
       to this function are blkEqn==23 and ptEqn==71, then the return value is
       2, which is the offset of ptEqn==71 into blkEqn==23.
   */
   int getBlkEqnOffset(int blkEqn, int ptEqn);

   /** Query whether specified point-equation is the first equation within
       a block-equation
   */
   bool isExactlyBlkEqn(int ptEqn);

   /** Return database of point-equations
    */
   std::map<int,int>* getPtEqns()
     {
       return( ptEqns_ );
     }

   /** Return database of block-equations
    */
   std::map<int,std::pair<int,int> >* getBlkEqns()
     {
       return( blkEqns_ );
     }

  private:
   /** copy constructor */
   PointBlockMap(const PointBlockMap& src);

   PointBlockMap& operator=(const PointBlockMap& src);

   std::map<int,int>* ptEqns_; //ptEqns_ maps point-equations to block-equations

   std::map<int,std::pair<int,int> >* blkEqns_;
   //blkEqns_ maps block-equations to point-equations and block-sizes
   //(the point-equation stored is the first or smallest point-equation that is
   //associated with the block-equation)

   int maxSize_;

   bool ptEqualBlk_; //if this is true, then all blkSizes are 1
  };//class PointBlockMap
}//namespace snl_fei
#endif // _snl_fei_PointBlockMap_hpp_

