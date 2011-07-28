/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Constraint_hpp_
#define _snl_fei_Constraint_hpp_

#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <fei_VectorSpace.hpp>
#include <snl_fei_RecordCollection.hpp>

#include <vector>

namespace snl_fei {

  /** container for constraint attributes */
  template<class RecordType>
  class Constraint {
  public:
    /** constructor */
    Constraint(int id=0, bool isPenaltyConstr=false);

    /** constructor */
    Constraint(int id,
               int constraintIDType,
               bool isSlave,
               bool isPenaltyConstr,
               int numIDs,
               const int* idTypes,
               const int* IDs,
               const int* fieldIDs,
               int offsetOfSlave,
               int offsetIntoSlaveField,
               const double* weights,
               double rhsValue,
               fei::VectorSpace* vspace);

    /** destructor */
    virtual ~Constraint();

    /** get constraint identifier */
    int getConstraintID() const { return( constraintID_ ); }

    /** set constraint identifier. power-users only */
    void setConstraintID(int id) { constraintID_ = id; }

    /** get the identifier-type that the fei uses to reference constraints */ 
    int getIDType() const { return( idType_ ); }

    snl_fei::RecordCollection* getRecordCollection() { return recordCollection_; }

    /** set the identifier-type that the fei uses to reference constraints.
      power-users only, this is a dangerous function with side-effects */ 
    void setIDType(int idType) { idType_ = idType; }

    /** query whether constraint is a penalty constraint */
    bool isPenalty() const { return( isPenalty_ ); }

    /** set whether constraint is a penalty constraint. Another dangerous
     function for power-users. */
    void setIsPenalty(bool isPenaltyConstr) { isPenalty_ = isPenaltyConstr; }

    /** get equation-number of constraint. (only valid if lagrange-multiplier)
    */
    int getEqnNumber() const { return( eqnNumber_ ); }

    /** set equation-number of constraint. (only valid if lagrange-multiplier)
    */
    void setEqnNumber(int eqn) { eqnNumber_ = eqn; }

    /** get block-equation number of constraint. (only valid if
     lagrange-multiplier */
    int getBlkEqnNumber() const { return( blkEqnNumber_ ); }

    /** set block-equation number of constraint. (only valid if
     lagrange-multiplier */
    void setBlkEqnNumber(int blkEqn) { blkEqnNumber_ = blkEqn; }


    /** intended for fei-implementation use only */
    RecordType getSlave() { return( slave_ ); }

    /** intended for fei-implementation use only */
    void setSlave(const RecordType& slv) { slave_ = slv; }

    /** if slave constraint, return constrained field-identifier of
        slaved mesh-object */
    int getSlaveFieldID() const { return( slaveField_ ); }

    /** if slave constraint, set constrained field-identifier of
        slaved mesh-object */
    void setSlaveFieldID(int f) { slaveField_ = f; }

    /** get offset of slaved field-component */
    int getOffsetIntoSlaveField() const { return( offsetIntoSlaveField_ ); }

    /** set offset of slaved field-component */
    void setOffsetIntoSlaveField(int offset) { offsetIntoSlaveField_ = offset; }


    /** get master mesh-objects */
    std::vector<int>& getMasters() { return( masters_ ); }

    /** get identifier-types of master mesh-objects */
    std::vector<int>& getMasterIDTypes() { return( masterIDTypes_ ); }

    /** get record-collections for masters */
    std::vector<snl_fei::RecordCollection*>& getMasterRecordCollections() { return masterRecordCollections_; }

    /** get field-identifiers of master mesh-objects */
    std::vector<int>& getMasterFieldIDs() { return( masterFields_ ); }

    /** get weight-coefficients of master mesh-objects */
    std::vector<double>& getMasterWeights() { return( masterWeights_ ); }


    /** get right-hand-side value of constraint */
    double getRHSValue() const { return( rhsValue_ ); }

    /** set right-hand-side value of constraint */
    void setRHSValue(double rhs) { rhsValue_ = rhs; }
 
    /** operator!= */
    bool operator!=(const Constraint<RecordType>& rhs);

    /** query whether connectivity is the same as specified constraint */
    bool structurallySame(const Constraint<RecordType>& rhs);

  private:
    Constraint(const Constraint<RecordType>& src);
    Constraint<RecordType>& operator=(const Constraint<RecordType>& src)
      {
        return(*this);
      }

    int constraintID_;
    int idType_;
    snl_fei::RecordCollection* recordCollection_;
    bool isPenalty_;

    int eqnNumber_;
    int blkEqnNumber_;

    RecordType slave_;
    int slaveField_;
    int offsetIntoSlaveField_;

    std::vector<int> masters_;
    std::vector<int> masterIDTypes_;
    std::vector<snl_fei::RecordCollection*> masterRecordCollections_;
    std::vector<int> masterFields_;
    std::vector<double> masterWeights_;

    double rhsValue_;

  };//class Constraint
} //namespace snl_fei

#include <snl_fei_Constraint.hpp>

//----------------------------------------------------------------------------
template<class RecordType>
inline snl_fei::Constraint<RecordType>::Constraint(int id, bool isPenaltyConstr)
  : constraintID_(id),
    idType_(0),
    recordCollection_(NULL),
    isPenalty_(isPenaltyConstr),
    eqnNumber_(-1),
    blkEqnNumber_(-1),
    slave_(),
    slaveField_(0),
    offsetIntoSlaveField_(0),
    masters_(),
    masterIDTypes_(),
    masterRecordCollections_(),
    masterFields_(),
    masterWeights_(),
    rhsValue_(0.0)
{
}

//----------------------------------------------------------------------------
template<class RecordType>
inline snl_fei::Constraint<RecordType>::Constraint(int id,
                                            int constraintIDType,
                                            bool isSlave,
                                            bool isPenaltyConstr,
                                            int numIDs,
                                            const int* idTypes,
                                            const int* IDs,
                                            const int* fieldIDs,
                                            int offsetOfSlave,
                                            int offsetIntoSlaveField,
                                            const double* weights,
                                            double rhsValue,
                                            fei::VectorSpace* vspace)
  : constraintID_(id),
    idType_(constraintIDType),
    recordCollection_(NULL),
    isPenalty_(isPenaltyConstr),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slave_(),
    slaveField_(0),
    offsetIntoSlaveField_(offsetIntoSlaveField),
    masters_(),
    masterIDTypes_(),
    masterRecordCollections_(),
    masterFields_(),
    masterWeights_(),
    rhsValue_(rhsValue)
{
}

//----------------------------------------------------------------------------
namespace snl_fei {
template<>
inline snl_fei::Constraint<fei::Record<int>*>::Constraint(int id,
                                            int constraintIDType,
                                            bool isSlave,
                                            bool isPenaltyConstr,
                                            int numIDs,
                                            const int* idTypes,
                                            const int* IDs,
                                            const int* fieldIDs,
                                            int offsetOfSlave,
                                            int offsetIntoSlaveField,
                                            const double* weights,
                                            double rhsValue,
                                            fei::VectorSpace* vspace)
  : constraintID_(id),
    idType_(constraintIDType),
    recordCollection_(NULL),
    isPenalty_(isPenaltyConstr),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slave_(),
    slaveField_(0),
    offsetIntoSlaveField_(offsetIntoSlaveField),
    masters_(),
    masterIDTypes_(),
    masterRecordCollections_(),
    masterFields_(),
    masterWeights_(),
    rhsValue_(rhsValue)
{
  int weightsOffset = 0;
  snl_fei::RecordCollection* recordCollection = NULL;
  vspace->getRecordCollection(idType_, recordCollection);
  recordCollection_ = recordCollection;
  for(int i=0; i<numIDs; ++i) {
    vspace->getRecordCollection(idTypes[i],recordCollection);
    masterRecordCollections_.push_back(recordCollection);

    vspace->addDOFs(fieldIDs[i], idTypes[i], 1, &(IDs[i]));
    int rec_local_id = recordCollection->getLocalID(IDs[i]);
    fei::Record<int>* rec = recordCollection->getRecordWithLocalID(rec_local_id);
    
    unsigned fieldSize = vspace->getFieldSize(fieldIDs[i]);

    if (isSlave && i == offsetOfSlave) {
      rec->hasSlaveDof(true);
      setSlave(rec);
      setSlaveFieldID(fieldIDs[i]);
      setOffsetIntoSlaveField(offsetIntoSlaveField);
      weightsOffset += fieldSize;
    }
    else {
      getMasters().push_back(rec_local_id);
      getMasterIDTypes().push_back(idTypes[i]);
      getMasterFieldIDs().push_back(fieldIDs[i]);

      if (weights != NULL) {
        for(unsigned j=0; j<fieldSize; ++j) {
          masterWeights_.push_back(weights[weightsOffset++]);
        }
      }
    }
  }
}

}//namespace snl_fei

//----------------------------------------------------------------------------
template<class RecordType>
inline snl_fei::Constraint<RecordType>::Constraint(const Constraint<RecordType>& src)
  : constraintID_(-1),
    idType_(0),
    isPenalty_(false),
    eqnNumber_(-1),
    blkEqnNumber_(-1), 
    slave_(),
    slaveField_(0),
    offsetIntoSlaveField_(0),
    masters_(),
    masterIDTypes_(),
    masterRecordCollections_(),
    masterFields_(),
    masterWeights_(),
    rhsValue_(0.0)
{
}

//----------------------------------------------------------------------------
template<class RecordType>
inline snl_fei::Constraint<RecordType>::~Constraint()
{
}

//----------------------------------------------------------------------------
template<class RecordType>
inline bool snl_fei::Constraint<RecordType>::operator!=(const snl_fei::Constraint<RecordType>& rhs)
{
  if (constraintID_ != rhs.constraintID_ ||
      idType_ != rhs.idType_ ||
      isPenalty_ != rhs.isPenalty_ ||
      eqnNumber_ != rhs.eqnNumber_ ||
      blkEqnNumber_ != rhs.blkEqnNumber_ ||
      slaveField_ != rhs.slaveField_ ||
      offsetIntoSlaveField_ != rhs.offsetIntoSlaveField_ ||
      rhsValue_ != rhs.rhsValue_) {
    return( true );
  }

  if (masters_ != rhs.masters_) return(true);

  if (masterIDTypes_ != rhs.masterIDTypes_) return(true);

  if (masterFields_ != rhs.masterFields_) return(true);

  if (masterWeights_ != rhs.masterWeights_) return(true);

  return(false);
}

//----------------------------------------------------------------------------
template<class RecordType>
inline bool snl_fei::Constraint<RecordType>::structurallySame(const Constraint<RecordType>& rhs)
{
  if (constraintID_ != rhs.constraintID_ ||
      idType_ != rhs.idType_ ||
      isPenalty_ != rhs.isPenalty_ ||
      eqnNumber_ != rhs.eqnNumber_ ||
      blkEqnNumber_ != rhs.blkEqnNumber_ ||
      slaveField_ != rhs.slaveField_ ||
      offsetIntoSlaveField_ != rhs.offsetIntoSlaveField_) {
    return( false );
  }

  if (masters_ != rhs.masters_) return(false);

  if (masterIDTypes_ != rhs.masterIDTypes_) return(false);

  if (masterFields_ != rhs.masterFields_) return(false);

  return(true);
}

#endif // _snl_fei_Constraint_hpp_

