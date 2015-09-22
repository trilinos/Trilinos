#ifndef __Panzer_Response_Residual_hpp__
#define __Panzer_Response_Residual_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

#include "Panzer_ResponseBase.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {

/** This class contains the data for construction
  * of a residual response. Note that the default version
  * is essentially empty.
  */
template <typename EvalT>
class Response_Residual : public ResponseBase {
public:
  Response_Residual(const std::string & responseName,
                     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof) :
    ResponseBase(responseName) {}
  virtual ~Response_Residual() {}
  virtual void scatterResponse() {}
  virtual void initializeResponse() {}
};

/** This is the response object used for calculation of the 
  * residual. This class uses the LOF to construct a ghosted residual
  * object. A user can uses class members to construct a compatible
  * residual object and then set it as the residual for this response to
  * fill.
  */
template < >
class Response_Residual<panzer::Traits::Residual> : public ResponseBase {
private:
  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;

  Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > residual_;
  mutable Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > ghostedResidual_;
      // mutable because of lazy construction (as needed)

public:
  Response_Residual(const std::string & responseName,
                     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof) 
     : ResponseBase(responseName) 
     , linObjFactory_(lof) {}
  virtual ~Response_Residual() {}

  /** Access the ghosted residual object. Note that this method will not return null.
    * When called for the first time this will use the LOF to construct a ghosted residual.
    */
  Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > getGhostedResidual() const; 

  /** Access the residual. This method can return null, but will only return the residual
    * class set by setResidual.
    */
  Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > getResidual() const; 

  /** Set the residual to use. If set to null, the internal residual will be lost. This
    * is assumed to be correctly sized.
    */
  void setResidual(const Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > & res);

  /** Build a correctly sized residual vector. This is a conenience, it wraps the 
    * linear object factory.
    */
  Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > allocateResidualVector() const; 

  // Functions inherited from ResponseBase
  virtual void initializeResponse() {}
  virtual void scatterResponse() {}
};

/** This is the response object used for calculation of the 
  * Jacobian. This class uses the LOF to construct a ghosted Jacobian
  * object. A user can uses class members to construct a compatible
  * Jacobian object and then set it as the Jacobian for this response to
  * fill.
  */
template < >
class Response_Residual<panzer::Traits::Jacobian> : public ResponseBase {
private:
  Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;

  Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > jacobian_;
  mutable Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > ghostedJacobian_;
      // mutable because of lazy construction (as needed)

public:
  Response_Residual(const std::string & responseName,
                     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof) 
     : ResponseBase(responseName) 
     , linObjFactory_(lof) {}
  virtual ~Response_Residual() {}

  /** Access the ghosted Jacobian object. Note that this method will not return null.
    * When called for the first time this will use the LOF to construct a ghosted Jacobian.
    */
  Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > getGhostedJacobian() const; 

  /** Access the Jacobian. This method can return null, but will only return the Jacobian
    * class set by setJacobian.
    */
  Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > getJacobian() const; 

  /** Set the Jacobian to use. If set to null, the internal Jacobian will be lost. This
    * is assumed to be correctly sized.
    */
  void setJacobian(const Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > & res);

  /** Build a correctly sized Jacobian. This is a conenience, it wraps the 
    * linear object factory.
    */
  Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > allocateJacobian() const; 

  // Functions inherited from ResponseBase
  virtual void initializeResponse() {}
  virtual void scatterResponse() {}
};

}

#endif
