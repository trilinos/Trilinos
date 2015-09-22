#ifndef _CRSet_h_
#define _CRSet_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/** Very simply struct-like class to hold data that defines a
Constraint Relation Set. A "set" of constraint relations implies that there may
be more than one, but current FEI usage (as of FEI version 2.0) results
in a "Set" always containing only one constraint relation. But it was deemed
unnecessary to change the name of this class.

It should also be noted that this container can also be used to hold the
definition for a "Slave variable" which is very similar to a constraint
relation, except that one of the constrained degrees of freedom is designated
as the "slave" -- it is defined to be a linear-combination of other degrees of
freedom plus a right-hand-side value (often zero).
*/

class CRSet {
 public:
  /** Default constructor. Does nothing but initialize pointer members to NULL.
   */
  CRSet();

  /** Destructor. Deletes arrays that were allocated during the life of this
      class instance.
  */
  ~CRSet();

  /** Identifier for this constraint relation. */
  int crID_;

  /** Number of nodes involved in this constraint relation. */
  int numNodes_;

  /** The node identifier of the node on which the "slaved" degree of freedom
      resides.
  */
  GlobalID slaveNodeID_;

  /** The field identifier of the slaved degree of freedom. */
  int slaveFieldID_;

  /** The offset into the field of the slaved degree of freedom. (Fields may 
      consist of several scalar degrees of freedom.)
  */
  int slaveOffset_;

   GlobalID** nodeIDs_;

   int* fieldIDs_;

   double* weights_;
   double* values_;
   double* penValues_;

 private:
   void deleteMemory();
};

#endif

