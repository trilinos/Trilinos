#ifndef __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
#define __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"

#include "Panzer_GlobalEvaluationData.hpp"

namespace panzer {

/** This class encapsulates the needs of a gather operation to do a halo
  * exchange.
  */
class ReadOnlyVector_GlobalEvaluationData : public GlobalEvaluationData {
public:

  //! Virtual d
  virtual ~ReadOnlyVector_GlobalEvaluationData() {}

  //! Is this object initialized
  virtual bool isInitialized() const = 0;

  /** For this class, this method does the halo exchange for the
    * vector.
    */
  virtual void globalToGhost(int mem) = 0;

  /** For this class, this method does nothing.
    */
  virtual void ghostToGlobal(int mem) {}

  //! Set the unique vector
  virtual void setUniqueVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & uniqueVector) = 0;

  //! Get the unique vector
  virtual Teuchos::RCP<const Thyra::VectorBase<double> > getUniqueVector() const = 0;

  //! Get the ghosted vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const = 0;
};

}

#endif
