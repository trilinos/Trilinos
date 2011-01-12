#ifndef stk_encr_GenericFunction_hpp
#define stk_encr_GenericFunction_hpp

#include <iostream>
#include <stk_percept/function/MDArray.hpp>
#include <stk_percept/function/internal/Dimensions.hpp>

namespace stk
{
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
    
    std::ostream &operator<<(std::ostream& out,  GenericFunction& func);
    //class NodalOp : public GenericFunction {};

  }
}

#endif
