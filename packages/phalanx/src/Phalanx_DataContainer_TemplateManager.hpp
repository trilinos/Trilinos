#ifndef PHX_DATA_CONTAINER_TEMPLATE_MANAGER_HPP
#define PHX_DATA_CONTAINER_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"
#include "Phalanx_DataContainer.hpp"

#include "Sacado_mpl_vector.hpp"
#include "boost/mpl/at.hpp"
#include "boost/mpl/placeholders.hpp"
using namespace boost::mpl::placeholders;


namespace PHX {
  
  template<class EvalT, typename Traits>
  class DataContainer_TemplateManager :
    public PHX::TemplateManager<typename boost::mpl::at<typename Traits::EvalToDataMap, EvalT>::type,
				PHX::DataContainerBase<Traits>,
				PHX::DataContainer<_, Traits> > {
    
  public:
    
    DataContainer_TemplateManager() {}
    
    ~DataContainer_TemplateManager() {}
    
  };
  
}
  
#endif 
