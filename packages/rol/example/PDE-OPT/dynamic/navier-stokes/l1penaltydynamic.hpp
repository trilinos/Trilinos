// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef L1DYNOBJECTIVE_HPP
#define L1DYNOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "../../TOOLS/ltiobjective.hpp"
#include "../../TOOLS/meshreader.hpp"
#include "../../TOOLS/pdevector.hpp"


template <class Real>
class L1_Dyn_Objective : public ROL::Objective<Real> {
	using size_type = typename std::vector<Real>::size_type; 
private:
  Real theta_, beta_;//, T_;
	const size_type  Nt_; 
  const std::vector<ROL::TimeStamp<Real>> ts_;
	const ROL::Ptr<ROL::PartitionedVector<Real>> zl_; 
	const ROL::Ptr<ROL::PartitionedVector<Real>> zu_; 
	
	ROL::PartitionedVector<Real> &partition ( ROL::Vector<Real>& x ) const {
    return static_cast<ROL::PartitionedVector<Real>&>(x);
  }

  const ROL::PartitionedVector<Real> &partition ( const ROL::Vector<Real>& x ) const {
    return static_cast<const ROL::PartitionedVector<Real>&>(x);
  }
	
	ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }


 public:
  L1_Dyn_Objective(ROL::ParameterList                         &parlist,
                   const std::vector<ROL::TimeStamp<Real>>    &timeStamp,
									 const ROL::Ptr<ROL::PartitionedVector<Real>>         &zl,
									 const ROL::Ptr<ROL::PartitionedVector<Real>>         &zu
    ) : Nt_(timeStamp.size()), ts_(timeStamp), zl_(zl), zu_(zu)
	    {
      theta_  = parlist.sublist("Reduced Dynamic Objective").sublist("Time Discretization").get("Theta",    1.0);
	    beta_   = parlist.sublist("Problem").get("L1 Control Cost", 1.e-2);

	  }
  
    // value
    Real value(const ROL::Vector<Real> &z,
               Real &tol 
		){
    
		  const ROL::PartitionedVector<Real> &zp = partition(z); 
		  const Real one(1);
      Real dt(0), val(0);
			bool isinf(false); 

      for (size_type k = 0; k < Nt_; ++k){//dynamic obj
				ROL::Ptr<const std::vector<Real>> zlk = getConstParameter(*zl_->get(k)); 
				ROL::Ptr<const std::vector<Real>> zuk = getConstParameter(*zu_->get(k)); 
			  ROL::Ptr<const std::vector<Real>> zpk = getConstParameter(*zp.get(k)); // ROL vector of partition/kth timestep
				
				for (size_type i = 0; i < zpk->size(); ++i){
					if ((*zpk)[i] < (*zlk)[i] || (*zpk)[i]>(*zuk)[i]){
						val = ROL::ROL_INF<Real>();  
						isinf=true; 
						break; 
					}
				}
			}
      
      if (!isinf) {
		    for (size_type k = 1; k < Nt_; ++k){//dynamic obj
		    	
		      ROL::Ptr<const std::vector<Real>> zpn = getConstParameter(*zp.get(k)); // ROL vector of partition/kth timestep
		      ROL::Ptr<const std::vector<Real>> zpo = getConstParameter(*zp.get(k-1)); // ROL vector of partition/kth timestep
	        dt = ts_[k].t[1] - ts_[k].t[0];

          for (size_type i = 0; i < zpn->size(); ++i){
		        val += dt*((one - theta_)*std::abs((*zpo)[i]) + theta_*std::abs((*zpn)[i])); 
		      } // end i for
		    }// end k for
			}
      return beta_*val;
   } //end value

	// prox
	void prox(ROL::Vector<Real>       &Pz, 
			      const ROL::Vector<Real> &z,  
						Real                t, 
						Real               &tol		)
	  {

		//partitioned vectors
		const ROL::PartitionedVector<Real> &zp  = partition(z); 
    ROL::PartitionedVector<Real> &Pzp = partition(Pz); 

		//constants
		const Real one(1), zero(0);
		Real hk(0), hkplus(0);
		Real l1param(0); 

		for (size_type k = 0; k<Nt_; ++k) {//dynamic part; 0->Nt inclusive
      ROL::Ptr<const std::vector<Real>> zlk = getConstParameter(*zl_->get(k)); 
			ROL::Ptr<const std::vector<Real>> zuk = getConstParameter(*zu_->get(k)); 
			ROL::Ptr<const std::vector<Real>> zk  = getConstParameter(zp[k]); //kth timestep
			ROL::Ptr<std::vector<Real>> Pzk = getParameter(Pzp[k]); 

			// update prox parameter
			if (k == 0){
        hkplus  = ts_[k+1].t[1] - ts_[k+1].t[0]; 
				l1param = t*hkplus*beta_*(one - theta_); 
			} else if (k == (Nt_ - 1)) {
				hk      = ts_[k].t[1] - ts_[k].t[0]; 
        l1param = t*hk*beta_*theta_; 
			} else {
				hk      = ts_[k].t[1] - ts_[k].t[0]; 
				hkplus  = ts_[k+1].t[1] - ts_[k+1].t[0]; 
        l1param = t*beta_*(hkplus*(one - theta_) + theta_*hk); 
			}
						
			for (size_type i = 0; i < zk->size(); ++i){
				if ((*zk)[i] > l1param) { (*Pzk)[i] = (*zk)[i] - l1param; }
				else if ((*zk)[i] < -l1param) {(*Pzk)[i] = (*zk)[i] + l1param;}
				else { (*Pzk)[i] = zero;}

				(*Pzk)[i] = std::min((*zuk)[i], std::max((*zlk)[i], (*Pzk)[i])); 
			} // end i loop

		}//end k for

	}//end prox	 
};

#endif
