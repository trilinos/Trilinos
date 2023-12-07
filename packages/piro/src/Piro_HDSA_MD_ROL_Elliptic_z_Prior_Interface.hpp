#ifndef PIRO_HDSA_MD_ROL_ELLIPTIC_Z_PRIOR_INTERFACE_HPP
#define PIRO_HDSA_MD_ROL_ELLIPTIC_Z_PRIOR_INTERFACE_HPP

#include "HDSA_MD_Elliptic_z_Prior_Interface.hpp"
#include "HDSA_ROL_Vector.hpp"

namespace Piro
{

template <class RealT>
class HDSA_MD_ROL_Elliptic_z_Prior_Interface : public HDSA::MD_Elliptic_z_Prior_Interface<RealT> {

private:
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > EzOp_;
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invEzOp_;
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > massOp_;
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invMassOp_;
  const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > cholMassOp_;

public:

  HDSA_MD_ROL_Elliptic_z_Prior_Interface(RealT & alpha_z,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > EzOp, 
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invEzOp,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > massOp,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > invMassOp,
         const Teuchos::RCP<const Thyra::LinearOpBase<RealT> > cholMassOp): 
        HDSA::MD_Elliptic_z_Prior_Interface<RealT>(alpha_z),
        EzOp_(EzOp),
        invEzOp_(invEzOp),
        massOp_(massOp),
        invMassOp_(invMassOp),
        cholMassOp_(cholMassOp)
  {}
  
  virtual ~HDSA_MD_ROL_Elliptic_z_Prior_Interface() {}

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis (from the base class)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  void Apply_E_z_Inverse_Implem(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, bool transpose) const {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    if(Teuchos::is_null(invEzOp_)) {
      z_out_rol.rol_vec->set(*z_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_out_rol.rol_vec).getVector();
      auto trans = transpose ? Thyra::TRANS : Thyra::NOTRANS;
      Thyra::apply(*invEzOp_, trans, *in, out.ptr(), (RealT)1, (RealT)0);
    }  
  }

  virtual void Apply_E_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    Apply_E_z_Inverse_Implem(z_out, z_in, false);
  }

  virtual void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    Apply_E_z_Inverse_Implem(z_out, z_in, true);
  }

  void Apply_E_z_Implem(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, bool transpose) const {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    if(Teuchos::is_null(EzOp_)) {
      z_out_rol.rol_vec->set(*z_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_out_rol.rol_vec).getVector();
      auto trans = transpose ? Thyra::TRANS : Thyra::NOTRANS;
      Thyra::apply(*EzOp_, trans, *in, out.ptr(), (RealT)1, (RealT)0);
    }  
  }

    virtual void Apply_E_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    Apply_E_z_Implem(z_out, z_in, false);
  }

  virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    Apply_E_z_Implem(z_out, z_in, true);
  }

  virtual void Apply_M_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    if(Teuchos::is_null(massOp_)) {
      z_out_rol.rol_vec->set(*z_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_out_rol.rol_vec).getVector();
      Thyra::apply(*massOp_, Thyra::NOTRANS, *in, out.ptr(), (RealT)1, (RealT)0);
    } 
  }

  virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const override {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    if(Teuchos::is_null(invMassOp_)) {
      z_out_rol.rol_vec->set(*z_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_out_rol.rol_vec).getVector();
      Thyra::apply(*invMassOp_, Thyra::NOTRANS, *in, out.ptr(), (RealT)1, (RealT)0);
    } 
  }

  virtual void Apply_CholM_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
    HDSA::ROL_Vector<RealT>& z_out_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(z_out);
    const HDSA::ROL_Vector<RealT>& z_in_rol = dynamic_cast<const HDSA::ROL_Vector<RealT>&>(z_in);
    if(Teuchos::is_null(cholMassOp_)) {
      z_out_rol.rol_vec->set(*z_in_rol.rol_vec);
    } else {
      Teuchos::RCP<const Thyra::VectorBase<RealT> > in = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_in_rol.rol_vec).getVector();
      Teuchos::RCP<Thyra::VectorBase<RealT> > out = dynamic_cast<ROL::ThyraVector<RealT> &>(*z_out_rol.rol_vec).getVector();
      Thyra::apply(*cholMassOp_, Thyra::TRANS, *in, out.ptr(), (RealT)1, (RealT)0);
    } 
  }

  void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> & samples) const override {
    int num_samples = samples.Number_of_Vectors();
    for(int i = 0; i < num_samples; i++)
    {
      auto sample = samples[i]->Clone();

      HDSA::ROL_Vector<RealT>& sample_rol = dynamic_cast<HDSA::ROL_Vector<RealT>&>(*sample);
      sample_rol.Randomize_Standard_Normal();

      auto chol_sample = samples[i]->Clone();
      Apply_CholM_z_Transpose(*chol_sample,*sample);
      Apply_E_z_Inverse(*sample,*chol_sample);
    } 
  }

};

}

#endif


