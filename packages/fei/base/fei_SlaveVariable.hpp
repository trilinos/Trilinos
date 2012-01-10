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

#ifndef _SlaveVariable_hpp_
#define _SlaveVariable_hpp_


/** Simple class to hold the information necessary to specify a slave variable
    in terms of a nodeID/fieldID/offsetIntoField and a list of master nodes with
    fields and coefficient-weights.
*/
class SlaveVariable {
 public:
  /** Default constructor */
  SlaveVariable()
    : nodeID_(-1), fieldID_(-1), offset_(0){
    masterNodes_ = new std::vector<GlobalID>; masterFields_ = new std::vector<int>;
    weights_ = new std::vector<double>;
  }

  /** Destructor */
  ~SlaveVariable() {delete masterNodes_; delete masterFields_; delete weights_;}

  GlobalID getNodeID() {return(nodeID_);}
  void setNodeID(GlobalID nid) {nodeID_ = nid;}

  int getFieldID() {return(fieldID_);}
  void setFieldID(int fid) {fieldID_ = fid;}

  int getFieldOffset() {return(offset_);}
  void setFieldOffset(int foff) {offset_ = foff;}

  const std::vector<GlobalID>* getMasterNodeIDs() {return(masterNodes_);}
  const std::vector<int>* getMasterFields() {return(masterFields_);}
  const std::vector<double>* getWeights() {return(weights_);}

  void addMasterNodeID(GlobalID masterNode) 
    {masterNodes_->push_back(masterNode);}

  void addMasterField(int masterField)
    {masterFields_->push_back(masterField);}

  void addWeight(double weight)
    {weights_->push_back(weight);}

 private:
  GlobalID nodeID_;
  int fieldID_;
  int offset_;

  std::vector<GlobalID>* masterNodes_;
  std::vector<int>* masterFields_;
  std::vector<double>* weights_;
};

#endif
