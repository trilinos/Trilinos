#ifndef _TSF_RDP_OPERATOR_H_
#define _TSF_RDP_OPERATOR_H_

//! TSF_RDP_Operator:  The Trilinos Real Double Precision Operator Class.
/*! The TSF_RDP_Operator class is a pure virtual class that specifies
    the required interfaces that any Trilinos-compliant real double precision 
    operator must implement.
*/
#include "TSF_Parameter.h"
#include "TSF_RDP_MultiVector.h"

class TSF_RDP_Operator {
    
  public:
  //! TSF_RDP_Operator Destructor.
  virtual ~TSF_RDP_Operator(void){};


  //! Allows user to pass list of parameters that can modify the behavior of the solver. 
  /*! The user can pass in one or more TSF_Parameter objects that encapsulate information
      for the solver to decode and either accept or reject as a valid parameter for the solver.
    \param In 
           NumParameters - Length of Parameters array.
    \param In 
           Parameters - An array of TSF_Parameter objects that will be scanned by the operator to determine
	   if the operator's behavior should be modified.
    \param Out 
           ParameterAccepted - An array of bools such that the ith element is true if the ith parameter in
	   Parameters was recognized by the operator.  Allows for tracking of unintended parameter behavior.

    \return Integer error code, set to 0 if successful.
           
  */
  virtual int SetParameters(int NumParameters, TSF_Parameter * Parameters, bool * ParameterAccepted) = 0;

  //! Performs any required setup that is needed by Apply(), must be called after SetParameters().
  /*!
    \return Integer error code, set to 0 if successful.
  */
  virtual int OperatorSetup(void) = 0;

  //! Apply operator to a given TSF_RDP_Multivector X, put results in Y.
  /*! Assuming OperatorSetup() is called, applies the operator.
    \param In 
           X - User supplied input multivector.
    \param Out 
           Y - Output multivector.

    \return Integer error code, set to 0 if successful.

  */
  virtual int Apply(TSF_RDP_MultiVector & X, TSF_RDP_MultiVector & Y) = 0;


};

#endif /* _TSF_RDP_OPERATOR_H_ */
