// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
