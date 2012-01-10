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

#ifndef _CRSet_h_
#define _CRSet_h_


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

