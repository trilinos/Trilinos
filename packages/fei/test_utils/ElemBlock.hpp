#ifndef _ElemBlock_h_
#define _ElemBlock_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

class ElemBlock {
 public:
   ElemBlock();
   ~ElemBlock();

   GlobalID blockID_;
   int numElements_;
   int numNodesPerElement_;
   int* numFieldsPerNode_;
   int** nodalFieldIDs_;
   GlobalID* elemIDs_;
   GlobalID** elemConn_;
   int numStiffRows_;
   int elemFormat_;
   double*** elemStiff_;
   double** elemLoad_;
   int numElemDOF_;
   int* elemDOFFieldIDs_;
   int interleaveStrategy_;
   int lumpingStrategy_;

 private:
   void deleteMemory();
};

#endif

