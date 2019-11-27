// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmoother2_hpp
#define ReferenceMeshSmoother2_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmoother1.hpp>

  namespace percept {

    /// A Jacobian based optimization smoother, eg. 1/A - 1/W (A = local current Jacobian, W is for original mesh)
    /// LOCAL patch version
    class ReferenceMeshSmoother2 : public ReferenceMeshSmoother1 {
     
    public:  
        
      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length 
      ReferenceMeshSmoother2(double max_edge_length_factor=1.0,
                                        int inner_iterations = 100,
                                        double cpu_time = 0.0, 
                                        double grad_norm =1.e-8,
                                        int parallel_iterations = 20)
        : ReferenceMeshSmoother1(max_edge_length_factor, inner_iterations, cpu_time, grad_norm, parallel_iterations)
      {}


    protected:

      virtual double run_one_iteration();


    };


  }

#endif
#endif

