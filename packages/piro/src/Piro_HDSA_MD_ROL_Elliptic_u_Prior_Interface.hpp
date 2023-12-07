#ifndef PIRO_HDSA_MD_ROL_ELLIPTIC_U_PRIOR_INTERFACE_HPP
#define PIRO_HDSA_MD_ROL_ELLIPTIC_U_PRIOR_INTERFACE_HPP

#include "HDSA_MD_Elliptic_u_Prior_Interface.hpp"
#include "HDSA_ROL_Vector.hpp"

namespace Piro
{

template <class RealT>
class HDSA_MD_ROL_Elliptic_u_Prior_Interface : public HDSA::MD_Elliptic_u_Prior_Interface<RealT> {

private:
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invEllOp_;
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > massOp_;

public:

  HDSA_MD_ROL_Elliptic_u_Prior_Interface(RealT & alpha_u, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invEllOp,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > massOp):
          HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u, random_number_generator),
        invEllOp_(invEllOp),
        massOp_(massOp)
  { }
  
  virtual ~HDSA_MD_ROL_Elliptic_u_Prior_Interface() {}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis (from the base class)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  void Apply_E_u_Inverse_Transpose_Implem(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, bool transpose) const {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    if(Teuchos::is_null(invEllOp_)) {
      u_out_rol.rol_vec->set(*u_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*u_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*u_out_rol.rol_vec).getVector();
      auto trans = transpose ? Thyra::TRANS : Thyra::NOTRANS;
      Thyra::apply(*invEllOp_, trans, *in, out.ptr(), (RealT)1, (RealT)0);
    }  
  }

  virtual void Apply_E_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const override {
    Apply_E_u_Inverse_Transpose_Implem(u_out, u_in, false);
  }

  virtual void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const override {
    Apply_E_u_Inverse_Transpose_Implem(u_out, u_in, true);
  }

  virtual void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const override {
    HDSA::ROL_Vector<RealT>& u_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(u_out);
    const HDSA::ROL_Vector<RealT>& u_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(u_in);
    if(Teuchos::is_null(massOp_)) {
      u_out_rol.rol_vec->set(*u_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*u_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*u_out_rol.rol_vec).getVector();
      Thyra::apply(*massOp_, Thyra::NOTRANS, *in, out.ptr(), (RealT)1, (RealT)0);
    } 
  }

};

}

#endif


