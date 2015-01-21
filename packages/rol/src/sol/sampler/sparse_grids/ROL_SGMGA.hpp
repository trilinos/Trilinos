// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_SGMGA_HPP
#define ROL_SGMGA_HPP

# include "ROL_SandiaRules.hpp"
# include "Teuchos_RCP.hpp"

# include <string>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include <cmath>

namespace ROL {

class sgmga {
private:
  Teuchos::RCP<SandiaRules> webbur;

public:

  sgmga(void){webbur = Teuchos::rcp(new SandiaRules());}

  double *sgmga_aniso_balance ( double alpha_max, 
				int dim_num, 
				double level_weight[] );

  void sgmga_aniso_normalize ( int option, 
			       int dim_num, 
			       double level_weight[] );

  void sgmga_importance_to_aniso ( int dim_num, 
				   double importance[], 
				   double level_weight[] );

  void sgmga_index ( int dim_num, 
		     double level_weight[], 
		     int level_max, 
		     int rule[], 
		     int point_num, 
		     int point_total_num, 
		     int sparse_unique_index[],
		     int growth[], 
		     int sparse_order[], 
		     int sparse_index[] );

  void sgmga_point ( int dim_num, 
		     double level_weight[], 
		     int level_max, 
		     int rule[], 
		     int np[], 
		     double p[], 
		     void ( *gw_compute_points[] ) ( int order, 
						     int np, 
						     double p[], 
						     double x[] ),
		     int point_num, 
		     int sparse_order[], 
		     int sparse_index[], 
		     int growth[], 
		     double sparse_point[] );

  void sgmga_product_weight ( int dim_num, 
			      int order_1d[], 
			      int order_nd, 
			      int rule[], 
			      int np[], 
			      double p[], 
			      void ( *gw_compute_weights[] ) ( int order, 
							       int np, 
							       double p[], 
							       double w[] ),
			      double weight_nd[] );

  int sgmga_size ( int dim_num, 
		   double level_weight[], 
		   int level_max, 
		   int rule[], 
		   int np[], 
		   double p[], 
		   void ( *gw_compute_points[] ) ( int order, 
						   int np, 
						   double p[], 
						   double x[] ),
		   double tol, 
		   int growth[] );

  int sgmga_size_total ( int dim_num, 
			 double level_weight[], 
			 int level_max, 
			 int rule[], 
			 int growth[] );

  void sgmga_unique_index ( int dim_num, 
			    double level_weight[], 
			    int level_max, 
			    int rule[], 
			    int np[], 
			    double p[], 
			    void ( *gw_compute_points[] ) ( int order, 
							    int np, 
							    double p[], 
							    double x[] ),
			    double tol, 
			    int point_num, 
			    int point_total_num, 
			    int growth[], 
			    int sparse_unique_index[] );

  void sgmga_vcn ( int n, 
		   double level_weight[], 
		   int x[], 
		   double q_min, 
		   double q_max, 
		   bool *more );

  double sgmga_vcn_coef ( int n, 
			  double level_weight[], 
			  int x[], 
			  double q_max );

  double sgmga_vcn_coef_naive ( int n, 
				double level_weight[], 
				int x_max[], 
				int x[], 
				double q_min, 
				double q_max );

  void sgmga_vcn_naive ( int n, 
			 double level_weight[], 
			 int x_max[], 
			 int x[], 
			 double q_min, 
			 double q_max, 
			 bool *more );

  void sgmga_vcn_ordered ( int dim_num, 
			   double level_weight[], 
			   int x_max[], 
			   int x[], 
			   double q_min, 
			   double q_max, 
			   bool *more );

  void sgmga_vcn_ordered_naive ( int dim_num, 
				 double level_weight[], 
				 int x_max[], 
				 int x[], 
				 double q_min, 
				 double q_max, 
				 bool *more );

  void sgmga_weight ( int dim_num, 
		      double level_weight[], 
		      int level_max, 
		      int rule[], 
		      int np[], 
		      double p[], 
		      void ( *gw_compute_weights[] ) ( int order, 
						       int np, 
						       double p[], 
						       double w[] ),
		      int point_num, 
		      int point_total_num, 
		      int sparse_unique_index[], 
		      int growth[], 
		      double sparse_weight[] );

  void sgmga_write ( int dim_num, 
		     double level_weight[], 
		     int rule[], 
		     int np[],
		     double p[], 
		     int point_num, 
		     double sparse_weight[], 
		     double sparse_point[],
		     std::string file_name );
};

} // namespace ROL

#include "ROL_SGMGADef.hpp"
#endif
