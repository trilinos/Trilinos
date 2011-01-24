#ifndef __Panzer_LinearObjFactory_hpp__
#define __Panzer_LinearObjFactory_hpp__

#include <map>

#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace panzer {

/** Abstract factory that builds the linear algebra 
  * objects required for the assembly.
  */
class LinearObjFactory {
public:

   /** get a vector meant to be filled by the phalanx assembly
     * this vector has ghosted entries and needs to be summed over
     * all the processors before it is actually used to solve with.
     */
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getGhostedVector() const = 0;

   /** Get a ghosted matrix.
     */
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getGhostedMatrix() const = 0;

   /** Get a globally distributed vector.
     */
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > getVector() const = 0;
   /** Get a globally distributed matrix.
     */
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > getMatrix() const = 0;

   /** Do the communication to fill a global matrix from a ghosted
     * matrix.
     */
   virtual void ghostToGlobalMatrix(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & ghostA, 
                                      const Teuchos::RCP<Thyra::LinearOpBase<double> > & A) const = 0;

   /** Do the communication to fill a global vector from a ghosted
     * matrix.
     */
   virtual void ghostToGlobalVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & ghostA, 
                                    const Teuchos::RCP<Thyra::MultiVectorBase<double> > & A) const = 0;

   /** Do the communication to fill a ghosted vector from a global
     * vector.
     */
   virtual void globalToGhostVector(const Teuchos::RCP<const Thyra::MultiVectorBase<double> > & A, 
                                    const Teuchos::RCP<Thyra::MultiVectorBase<double> > & ghostA) const = 0;
};

}

#endif
