// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_GenericFunction_hpp
#define stk_encr_GenericFunction_hpp

#include <iostream>
#include <percept/function/MDArray.hpp>
#include <percept/function/internal/Dimensions.hpp>

  namespace percept
  {

    class GenericFunction
    {
    public:
      GenericFunction(Dimensions domain_dimensions = Dimensions(),
                      Dimensions codomain_dimensions = Dimensions()) :
        m_domain_dimensions(domain_dimensions),
        m_codomain_dimensions(codomain_dimensions), m_spatialOperator(false)
      {
      }
      virtual ~GenericFunction() {}

      /** Evaluate the function on it's domain returning result in codomain.  
       */
      virtual void operator()(MDArray& domain, MDArray& codomain, double time = 0.0)=0;

      bool isSpatialOperator() { return m_spatialOperator; }
      void setIsSpatialOperator(bool so) { m_spatialOperator=so; }

      Dimensions getDomainDimensions() {return m_domain_dimensions; }
      Dimensions getCodomainDimensions() { return m_codomain_dimensions; }
      MDArray getNewDomain()   
      { 
        return Intrepid::FieldContainer<double>( Teuchos::Array<int>(m_domain_dimensions.begin(),   m_domain_dimensions.end()));   
      }

      MDArray getNewCodomain() 
      { 
        return Intrepid::FieldContainer<double>( Teuchos::Array<int>(m_codomain_dimensions.begin(), m_codomain_dimensions.end())); 
      }
      MDArray getNewCodomain() const
      { 
        return Intrepid::FieldContainer<double>( Teuchos::Array<int>(m_codomain_dimensions.begin(), m_codomain_dimensions.end())); 
      }

      static MDArray getNewMDArray(const Dimensions dims) 
      { 
        return Intrepid::FieldContainer<double>( Teuchos::Array<int>(dims.begin(),   dims.end()));   
      }
    protected:
      Dimensions m_domain_dimensions;  // size() gives rank, each entry gives dimension, e.g. {3,3} for a rank-2 3D tensor
      Dimensions m_codomain_dimensions;
      bool m_spatialOperator;
    };
    
#ifndef SWIG
    std::ostream &operator<<(std::ostream& out,  GenericFunction& func);
#endif
    //class NodalOp : public GenericFunction {};

  }

#endif
