#ifndef __Panzer_WriteVector_GlobalEvaluationData_hpp__
#define __Panzer_WriteVector_GlobalEvaluationData_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"

#include "Panzer_GlobalEvaluationData.hpp"

namespace panzer {

/** This class encapsulates the needs of a gather operation to do a halo
  * exchange.
  */
class WriteVector_GlobalEvaluationData : public GlobalEvaluationData {
public:

  //! when you gho from ghost to global, combine with a particular mode
  typedef enum { CM_Sum, CM_Max, CM_Min, CM_Insert } CombineMode;

  //! Default constructor, set combine mode to sum right away
  WriteVector_GlobalEvaluationData() : combineMode_(CM_Sum) {}

  //! Virtual d
  virtual ~WriteVector_GlobalEvaluationData() {}

  //! Allow the user to set the combine mode (at any time)
  void setCombineMode(CombineMode cm) { combineMode_ = cm; }

  //! Get the combine mode, to be used by sub classes
  CombineMode getCombineMode() const { return combineMode_; }

  //! Is this object initialized
  virtual bool isInitialized() const = 0;


  /** For this class, this method does nothing.
    */
  virtual void globalToGhost(int mem) { }

  /** For this class, this method does the halo exchange for the
    * vector.
    */
  virtual void ghostToGlobal(int mem) = 0;

  //! Set the owned vector
  virtual void setOwnedVector(const Teuchos::RCP<Thyra::VectorBase<double> > & ownedVector) = 0;

  //! Get the owned vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getOwnedVector() const = 0;

  //! Get the ghosted vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const = 0;

private: 

  CombineMode combineMode_;
};

}

#endif
