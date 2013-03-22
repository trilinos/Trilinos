#ifndef _CommNodeSet_h_
#define _CommNodeSet_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

class CommNodeSet {
 public:
   CommNodeSet();
   ~CommNodeSet();

   int numNodes_;
   GlobalID* nodeIDs_;
   int** procs_;
   int* procsPerNode_;
 private:
   void deleteMemory();
};

#endif

