#ifndef PIRO_HDSA_MD_ROL_DATA_INTERFACE_HPP
#define PIRO_HDSA_MD_ROL_DATA_INTERFACE_HPP

#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_ROL_Vector.hpp"

namespace Piro
{

template <class RealT>
class HDSA_MD_ROL_Data_Interface : public HDSA::MD_Data_Interface<RealT> {

private:
  const ROL::Ptr<ROL::Vector<RealT> > & u_rol_opt_;
  const ROL::Ptr<ROL::Vector<RealT> > & z_rol_opt_;
  mutable HDSA::Ptr<HDSA::MultiVector<RealT> > d_data_;
  mutable HDSA::Ptr<HDSA::MultiVector<RealT> > z_data_;


public:
  HDSA_MD_ROL_Data_Interface(const ROL::Ptr<ROL::Vector<RealT> > & u_opt, const ROL::Ptr<ROL::Vector<RealT> > & z_opt):
    u_rol_opt_(u_opt), 
    z_rol_opt_(z_opt)
    {  
      d_data_ = HDSA::makePtr<HDSA::MultiVector<RealT> >();
      z_data_ = HDSA::makePtr<HDSA::MultiVector<RealT> >();
    }

  virtual ~HDSA_MD_ROL_Data_Interface()
  { }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const {
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(u_rol_opt_);
    return u_opt;
  }
  
  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const {
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z_rol_opt_);
    return z_opt;
  }
  
  void Z_Data_push_back(const ROL::Ptr<ROL::Vector<RealT> > & z_data_rol) const {
    HDSA::Ptr<HDSA::Vector<RealT> > z_data = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(z_data_rol);
    z_data_->push_back(z_data);
  }

  void Y_Data_push_back(const ROL::Ptr<ROL::Vector<RealT> > & d_data_rol) const {
    HDSA::Ptr<HDSA::Vector<RealT> > d_data = HDSA::makePtr<HDSA::ROL_Vector<RealT> >(d_data_rol);
    d_data_->push_back(d_data);
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const override {
    return d_data_;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const override {
    return z_data_;
  }  
};

} //namespace Piro


#endif

