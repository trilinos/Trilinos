#ifndef _SlaveVariable_hpp_
#define _SlaveVariable_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
