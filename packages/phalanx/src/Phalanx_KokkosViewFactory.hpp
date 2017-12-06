#ifndef PHALANX_KOKKOS_VIEW_FACTORY_HPP
#define PHALANX_KOKKOS_VIEW_FACTORY_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_any.hpp"
#include "Sacado.hpp"
#include <vector>

namespace PHX {
  
  template<typename ScalarT, typename Device>
  class KokkosViewFactory {
  public:
    static PHX::any buildView(const PHX::FieldTag& t, const std::vector<PHX::index_size_type>& extended_dimensions = std::vector<PHX::index_size_type>(0));
  };

  // *********************************************
  // Default implementation.  NOTE: for DFad types we implement
  // partial specializations below.
  // *********************************************
  template<typename ScalarT, typename Device>
  PHX::any 
  KokkosViewFactory<ScalarT,Device>::buildView(const PHX::FieldTag& t,
					       const std::vector<PHX::index_size_type>& )
  {
    PHX::any a;
    const PHX::DataLayout& dl = t.dataLayout();

    if (dl.rank() == 1)
      a = PHX::View<ScalarT*>(t.identifier(),
			      dl.dimension(0));
    else if (dl.rank() == 2)
      a = PHX::View<ScalarT**>(t.identifier(),
			       dl.dimension(0),
			       dl.dimension(1));
    else if (dl.rank() == 3)
      a = PHX::View<ScalarT***>(t.identifier(),
				dl.dimension(0),
				dl.dimension(1),
				dl.dimension(2));
    else if (dl.rank() == 4)
      a = PHX::View<ScalarT****>(t.identifier(),
				 dl.dimension(0),
				 dl.dimension(1),
				 dl.dimension(2),
				 dl.dimension(3));
    else if (dl.rank() == 5)
      a = PHX::View<ScalarT*****>(t.identifier(),
				  dl.dimension(0),
				  dl.dimension(1),
				  dl.dimension(2),
				  dl.dimension(3),
				  dl.dimension(4));
    else if (dl.rank() == 6)
      a = PHX::View<ScalarT******>(t.identifier(),
				   dl.dimension(0),
				   dl.dimension(1),
				   dl.dimension(2),
				   dl.dimension(3),
				   dl.dimension(4),
				   dl.dimension(5));
    else if (dl.rank() == 7)
      a = PHX::View<ScalarT*******>(t.identifier(),
				    dl.dimension(0),
				    dl.dimension(1),
				    dl.dimension(2),
				    dl.dimension(3),
				    dl.dimension(4),
				    dl.dimension(5),
				    dl.dimension(6));
    return a;
  }
  
  // *********************************************
  // Sacado::Fad::DFad Partial Specialization
  // *********************************************
  template<typename ScalarT, typename Device>
  class KokkosViewFactory<Sacado::Fad::DFad<ScalarT>,Device> {
  public:
    static 
    PHX::any 
    buildView(const PHX::FieldTag& t,
	      const std::vector<PHX::index_size_type>& derivative_dimensions)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(derivative_dimensions.size() != 1,std::runtime_error,
				 "Error in PHX::KokkosViewFactory::buildView() the Sacado Fad type requires a runtime specified array dimension size of the number of derivitives.  Please set this using the function \"setKokkosExtendedDataTypeDimensions()\" on the field manager!");

      PHX::any a;
      const PHX::DataLayout& dl = t.dataLayout();
      // DFad type contains a hidden dimension of the size of the number
      // of derivatives.  We add one to this for the value of the
      // function.
      const PHX::index_size_type hDim = derivative_dimensions[0] + 1;
      
      if (dl.rank() == 1)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>*>(t.identifier(),
						   dl.dimension(0),
						   hDim);
      else if (dl.rank() == 2)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>**>(t.identifier(),
						    dl.dimension(0),
						    dl.dimension(1),
						    hDim);
      else if (dl.rank() == 3)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>***>(t.identifier(),
						     dl.dimension(0),
						     dl.dimension(1),
						     dl.dimension(2),
						     hDim);
      else if (dl.rank() == 4)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>****>(t.identifier(),
						      dl.dimension(0),
						      dl.dimension(1),
						      dl.dimension(2),
						      dl.dimension(3),
						      hDim);
      else if (dl.rank() == 5)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>*****>(t.identifier(),
						       dl.dimension(0),
						       dl.dimension(1),
						       dl.dimension(2),
						       dl.dimension(3),
						       dl.dimension(4),
						       hDim);
      else if (dl.rank() == 6)
	a = PHX::View<Sacado::Fad::DFad<ScalarT>******>(t.identifier(),
							dl.dimension(0),
							dl.dimension(1),
							dl.dimension(2),
							dl.dimension(3),
							dl.dimension(4),
							dl.dimension(5),
							hDim);
      else if (dl.rank() == 7)
      	a = PHX::View<Sacado::Fad::DFad<ScalarT>*******>(t.identifier(),
							 dl.dimension(0),
							 dl.dimension(1),
							 dl.dimension(2),
							 dl.dimension(3),
							 dl.dimension(4),
							 dl.dimension(5),
							 dl.dimension(6),
							 hDim);
      
      
      return a;
    }
  };
  
  // *********************************************
  // Sacado::ELRCacheFad::DFad Partial Specialization
  // *********************************************
  template<typename ScalarT, typename Device>
  class KokkosViewFactory<Sacado::ELRCacheFad::DFad<ScalarT>,Device> {
  public:
    static 
    PHX::any 
    buildView(const PHX::FieldTag& t,
	      const std::vector<PHX::index_size_type>& derivative_dimensions)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(derivative_dimensions.size() != 1,std::runtime_error,
				 "Error in PHX::KokkosViewFactory::buildView() the Sacado Fad type requires a runtime specified array dimension size of the number of derivitives.  Please set this using the function \"setKokkosExtendedDataTypeDimensions()\" on the field manager!");

      PHX::any a;
      const PHX::DataLayout& dl = t.dataLayout();
      // DFad type contains a hidden dimension of the size of the number
      // of derivatives.  We add one to this for the value of the
      // function.
      const PHX::index_size_type hDim = derivative_dimensions[0] + 1;
      
      if (dl.rank() == 1)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>*>(t.identifier(),
							   dl.dimension(0),
							   hDim);
      else if (dl.rank() == 2)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>**>(t.identifier(),
							    dl.dimension(0),
							    dl.dimension(1),
							    hDim);
      else if (dl.rank() == 3)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>***>(t.identifier(),
							     dl.dimension(0),
							     dl.dimension(1),
							     dl.dimension(2),
							     hDim);
      else if (dl.rank() == 4)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>****>(t.identifier(),
							      dl.dimension(0),
							      dl.dimension(1),
							      dl.dimension(2),
							      dl.dimension(3),
							      hDim);
      else if (dl.rank() == 5)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>*****>(t.identifier(),
							       dl.dimension(0),
							       dl.dimension(1),
							       dl.dimension(2),
							       dl.dimension(3),
							       dl.dimension(4),
							       hDim);
      else if (dl.rank() == 6)
	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>******>(t.identifier(),
								dl.dimension(0),
								dl.dimension(1),
								dl.dimension(2),
								dl.dimension(3),
								dl.dimension(4),
								dl.dimension(5),
								hDim);
      else if (dl.rank() == 7)
      	a = PHX::View<Sacado::ELRCacheFad::DFad<ScalarT>*******>(t.identifier(),
								 dl.dimension(0),
								 dl.dimension(1),
								 dl.dimension(2),
								 dl.dimension(3),
								 dl.dimension(4),
								 dl.dimension(5),
								 dl.dimension(6),
								 hDim);
      
      
      return a;
    }
  };

  
  // *********************************************
  // Sacado::Fad::SLFad Partial Specialization
  // *********************************************
  template<typename ScalarT, typename Device, int N>
  class KokkosViewFactory<Sacado::Fad::SLFad<ScalarT,N>,Device> {
  public:
    static 
    PHX::any 
    buildView(const PHX::FieldTag& t,
	      const std::vector<PHX::index_size_type>& derivative_dimensions)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(derivative_dimensions.size() != 1,std::runtime_error,
				 "Error in PHX::KokkosViewFactory::buildView() the Sacado Fad type requires a runtime specified array dimension size of the number of derivitives.  Please set this using the function \"setKokkosExtendedDataTypeDimensions()\" on the field manager!");

      PHX::any a;
      const PHX::DataLayout& dl = t.dataLayout();
      // DFad type contains a hidden dimension of the size of the number
      // of derivatives.  We add one to this for the value of the
      // function.
      const PHX::index_size_type hDim = derivative_dimensions[0] + 1;

      TEUCHOS_TEST_FOR_EXCEPTION((derivative_dimensions[0] > N), std::logic_error,
				 "Error: The MDField \"" << t.identifier()
				 << "\" has an SLFAD compile time derivative size limit of " 
				 << N << ", but the runtime requested derivative size "
				 << derivative_dimensions[0] << ", is larger.");
      
      if (dl.rank() == 1)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>*>(t.identifier(),
						      dl.dimension(0),
						      hDim);
      else if (dl.rank() == 2)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>**>(t.identifier(),
						       dl.dimension(0),
						       dl.dimension(1),
						       hDim);
      else if (dl.rank() == 3)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>***>(t.identifier(),
							dl.dimension(0),
							dl.dimension(1),
							dl.dimension(2),
							hDim);
      else if (dl.rank() == 4)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>****>(t.identifier(),
							 dl.dimension(0),
							 dl.dimension(1),
							 dl.dimension(2),
							 dl.dimension(3),
							 hDim);
      else if (dl.rank() == 5)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>*****>(t.identifier(),
							  dl.dimension(0),
							  dl.dimension(1),
							  dl.dimension(2),
							  dl.dimension(3),
							  dl.dimension(4),
							  hDim);
      else if (dl.rank() == 6)
	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>******>(t.identifier(),
							   dl.dimension(0),
							   dl.dimension(1),
							   dl.dimension(2),
							   dl.dimension(3),
							   dl.dimension(4),
							   dl.dimension(5),
							   hDim);
      else if (dl.rank() == 7)
      	a = PHX::View<Sacado::Fad::SLFad<ScalarT,N>*******>(t.identifier(),
							    dl.dimension(0),
							    dl.dimension(1),
							    dl.dimension(2),
							    dl.dimension(3),
							    dl.dimension(4),
							    dl.dimension(5),
							    dl.dimension(6),
							    hDim);
      
      
      return a;
    }
  };

}

#endif
