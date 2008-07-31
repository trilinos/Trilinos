#ifndef PHX_TRAITS_HPP
#define PHX_TRAITS_HPP

// mpl (Meta Programming Library) templates
#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"
#include "boost/mpl/map.hpp"
#include "boost/mpl/find.hpp"

// traits Base Class
#include "Phalanx_Traits_Base.hpp"

// Include User Data Types
#include "Phalanx_ConfigDefs.hpp" // for std::vector
#include "AlgebraicTypes.hpp"
#include "CellData.hpp"
#include "Phalanx_Allocator_New.hpp"

namespace PHX {

  /*! \brief Traits class used for testing.
    
  */
  struct MyTraits : public PHX::TraitsBase {
    
    // ******************************************************************
    // *** Scalar Types
    // ******************************************************************
    
    // Scalar types we plan to use
    typedef double RealType;
    typedef Sacado::Fad::DFad<double> FadType;
    
    // ******************************************************************
    // *** Evaluation Types
    // ******************************************************************
    struct Residual { typedef RealType ScalarT; };
    struct Jacobian { typedef FadType ScalarT;  };
    typedef Sacado::mpl::vector<Residual, Jacobian> EvalTypes;

    // ******************************************************************
    // *** Data Types
    // ******************************************************************
    
    // Create the data types for each scalar type
    
    // <double>
    typedef Sacado::mpl::vector< RealType, 
				 MyVector<RealType>,
				 MyTensor<RealType> 
    > ResidualDataTypes;
  
    // Fad<double, double>
    typedef Sacado::mpl::vector< FadType,
				 MyVector<FadType>,
				 MyTensor<FadType> 
    > JacobianDataTypes;

    // Maps the key EvalType a vector of DataTypes
    typedef boost::mpl::map<
      boost::mpl::pair<Residual, ResidualDataTypes>,
      boost::mpl::pair<Jacobian, JacobianDataTypes>
    >::type EvalToDataMap;

    // ******************************************************************
    // *** Algebraic Types
    // ******************************************************************

    struct MY_SCALAR { static const std::string name; };
    struct MY_VECTOR { static const std::string name; };
    struct MY_TENSOR { static const std::string name; };

    // ******************************************************************
    // *** Connectors Types
    // ******************************************************************

    // Maps the key DataType to the value AlgebraicTypes
    typedef boost::mpl::map<
      boost::mpl::pair< RealType          , MY_SCALAR>,
      boost::mpl::pair< MyVector<RealType>, MY_VECTOR>,
      boost::mpl::pair< MyTensor<RealType>, MY_TENSOR>,
      boost::mpl::pair< FadType           , MY_SCALAR>,
      boost::mpl::pair< MyVector<FadType> , MY_VECTOR>,
      boost::mpl::pair< MyTensor<FadType> , MY_TENSOR>
    >::type DataToAlgebraicMap;
    
    // ******************************************************************
    // *** Allocator Type
    // ******************************************************************
    typedef PHX::NewAllocator Allocator;

    // ******************************************************************
    // *** User Defined Object Passed in for Evaluation Method
    // ******************************************************************
    typedef std::vector<CellData>& EvalData;
    typedef void* PreEvalData;
    typedef void* PostEvalData;

  };
 
  // ******************************************************************
  // ******************************************************************
  // Specializations of the DataTypeInfo found in 
  // PHX::TraitsBase::DataTypeInfo
  // ******************************************************************
  // ******************************************************************
  /*
  template<>
  struct MyTraits::DataTypeInfo<MyTraits::RealType>
  {
  typedef MyTraits::RealType scalar_type;
  typedef MyTraits::RealType algebric_type;
  };
  */

  // ******************************************************************
  // ******************************************************************
  // Debug strings.
  // 1. Initialize the name member for Algebric Types
  // 2. Specialize the Scalar types and Data types for the TypeString
  //    object in the PHX::TraitsBase class.
  // ******************************************************************
  // ******************************************************************

  // Define the string names for the Algebraic Types
  const std::string MyTraits::MY_SCALAR::name = "Scalar";
  const std::string MyTraits::MY_VECTOR::name = "Vector";
  const std::string MyTraits::MY_TENSOR::name = "Tensor";

  // Evaluation Types
  template<>
  struct MyTraits::TypeString<MyTraits::Residual> 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTraits::Residual>::value = 
    "Residual";

  template<>
  struct MyTraits::TypeString<MyTraits::Jacobian> 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTraits::Jacobian>::value = 
    "Jacobian";

  // Data Types
  template<>
  struct MyTraits::TypeString<MyTraits::RealType> 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTraits::RealType>::value = 
    "double";

  template<>
  struct MyTraits::TypeString<MyTraits::FadType> 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTraits::FadType>::value = 
    "Sacado::Fad::DFad<double>";

  template<>
  struct MyTraits::TypeString< MyVector<MyTraits::RealType> > 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyVector<MyTraits::RealType> >::value = 
    "MyVector<double>";

  template<>
  struct MyTraits::TypeString< MyVector<MyTraits::FadType> > 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyVector<MyTraits::FadType> >::value = 
    "MyVector< Sacado::Fad::DFad<double> >";

  template<>
  struct MyTraits::TypeString<MyTensor<MyTraits::RealType> > 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTensor<MyTraits::RealType> >::value = 
    "MyTensor<double>";

  template<>
  struct MyTraits::TypeString< MyTensor<MyTraits::FadType> > 
  { static const std::string value; };
  const std::string MyTraits::TypeString<MyTensor<MyTraits::FadType> >::value = 
    "MyTensor< Sacado::Fad::DFad<double> >";

}

#endif
