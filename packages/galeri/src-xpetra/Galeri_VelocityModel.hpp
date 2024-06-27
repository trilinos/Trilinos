// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_VELOCITYMODEL_HPP
#define GALERI_VELOCITYMODEL_HPP

// velocity models for Helmholtz

namespace Galeri {

  namespace Xpetra {

    template<typename Scalar,typename LocalOrdinal,typename GlobalOrdinal>
    class VelocityModel {

    public:

      VelocityModel(LocalOrdinal dimension,LocalOrdinal model) {
	dim_=dimension;
	modelType_=model;
      }

      VelocityModel() {
	dim_=3;
	modelType_=0;
      }

      void setDim(LocalOrdinal dimension) {
	dim_=dimension;
      }

      void setModel(LocalOrdinal model) {
	modelType_=model;
      }

      Scalar getVelocity(double x,double y,double z);

    private:

      typedef Scalar        SC;
      typedef LocalOrdinal  LO;

      // dimension and model type
      LO dim_, modelType_;

    };

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal>
    Scalar VelocityModel<Scalar,LocalOrdinal,GlobalOrdinal>::getVelocity(double x, double y, double z) {
      
      Scalar c=(Scalar)1.0;
      if(dim_==2) {
	if(modelType_==0) {
	  // constant wavespeed
	  c=(Scalar)1.0;
	}
	else if(modelType_==1) {
	  // converging lens
	  c=(Scalar)(4.0/3.0)*(1.0-0.5*exp(-32.0*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))));
	}
	else if(modelType_==2) {
	  // wedge
	  c=(Scalar)1.0;
	  if(x<0.8-0.2*y)
	    c=(Scalar)0.7;
	  if(x<0.4+0.1*y)
	    c=(Scalar)1.3;
	}
      }
      else { // if dim_==3
	if(modelType_==0) {
	  // constant wavespeed
	  c=(Scalar)1.0;
	}
	else if(modelType_==1) {
	  // converging lens
	  c=(Scalar)(4.0/3.0)*(1.0-0.5*exp(-32.0*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))));
	}
	else if(modelType_==2) {
	  // wedge
	  c=(Scalar)1.0;
	  if(x<0.8-0.2*z)
	    c=(Scalar)0.7;
	  if(x<0.4+0.1*z)
	    c=(Scalar)1.3;
	}
      }
      return c;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_HELMHOLTZ3DPROBLEM_HPP
