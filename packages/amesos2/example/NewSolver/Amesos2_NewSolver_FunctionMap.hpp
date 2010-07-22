/**
   \file   NewSolver_FunctionMap.hpp
   \author John Doe <jd@sandia.gov>
   \date   Fri Jul  9 14:48:16 CDT 2010
   
   \brief  Template for providing a mechanism to map function calls to the
           correct Solver function based on the scalar type of Matrices and
           MultiVectors
*/

#ifndef AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
#define AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP

#include <complex>

#include "Amesos2_FunctionMap.hpp"
#include "Amesos2_MatrixAdapter.hpp"
#include "NewSolver_TypeMap.hpp"


/*
 * External definitions of the NewSolver functions in namespace New_Solver
 */
namespace New_Solver {
extern "C" {

/*
 * Following is an example class.  Actual class would be #include'd from the
 * solver package library
 */
class Matrix;

/*
 * Include solver headers for different data types in different namespaces.
 */
namespace S {
// include here single-precision real definition headers
}

namespace D {
// include here double-precision real definition headers
}

namespace C {
// include here single-precision complex definition headers
}

namespace Z {
// include here double-precision complex definition headers
}

} // end extern "C"
} // end namespace New_Solver


namespace Amesos {

/** 
 * Helper class which passes on function calls to the appropriate NewSolver
 * function based on the type of its scalar template argument.
 *
 * Some solver libraries have solver and matrix builder functions defined
 * based on data type.  One function for complex, one for double precision
 * complex, another for \c float , and yet another for \c double.  To work
 * elegantly with the Amesos::NewSolver interface we want to be able to
 * perform a single function call which is appropriate for the scalar type of
 * the Matrix and MultiVectors that we are working with.  The \c FunctionMap
 * class provides that capability.
 *
 * The class template is specialized for each data type that NewSolver supports,
 * and errors are thrown for other data types.
 */
template <typename Scalar>
struct FunctionMap<NewSolver,Scalar>
{

  /**
   * \brief Binds to the appropriate NewSolver solver driver based on data type.
   * 
   * \throw std::runtime_error If no specialization of this type exists for a
   *        particular scalar type
   */
  static void solve(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }


  /**
   * \brief Binds to the appropriate NewSolver symbolic factorization method
   * based on data type.
   * 
   * \throw std::runtime_error If no specialization of this type exists for a
   *        particular scalar type
   */
  static void symbolic_fact(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }


  /**
   * \brief Binds to the appropriate NewSolver numeric factorization method
   * based on data type.
   * 
   * \throw std::runtime_error If no specialization of this type exists for a
   *        particular scalar type
   */
  static void numeric_fact(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }
  

  /** \brief Creates a NewSolver CCS matrix using the appropriate function
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompCol_Matrix(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }


  /** \brief Creates a NewSolver CRS matrix using the appropriate function
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_CompRow_Matrix(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }


  /** \brief Creates a NewSolver Dense Matrix using the appropriate NewSolver
   *         function.
   * 
   * \throw std::runtime_error If there is no specialization of this type for
   *        the Scalar type
   */
  static void create_Dense_Matrix(
    /* args */
    )
    {
      TEST_FOR_EXCEPTION( true,
        std::runtime_error,
        "NewSolver does not support the data type");
    }
};


/* ==================== Specializations ====================
 *
 * Provide a specialization for the template class which binds each function
 * to the appropriate NewSolver package function call for the given data
 * type.
 * 
 */
template <>
struct FunctionMap<NewSolver,float>
{
  static void solve( /* args */ )
    {
      New_Solver::S::single_solve( /* args */ );
    }

  static void symbolic_fact( /* args */ )
    {
      New_Solver::S::single_symbolic_fact( /* args */ );
    }

  static void numeric_fact( /* args */ )
    {
      New_Solver::S::single_numeric_fact( /* args */ );
    }

  static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::S::single_create_CompCol_Matrix( /* args */ );
    }

  static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::S::single_create_CompRow_Matrix( /* args */ );
    }

  static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::S::single_create_Dense_Matrix( /* args */ );
    }
};


template <>
struct FunctionMap<NewSolver,double>
{
  static void solve( /* args */ )
    {
      New_Solver::D::double_solve( /* args */ );
    }

  static void symbolic_fact( /* args */ )
    {
      New_Solver::D::double_symbolic_fact( /* args */ );
    }

  static void numeric_fact( /* args */ )
    {
      New_Solver::D::double_numeric_fact( /* args */ );
    }

  static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::D::double_create_CompCol_Matrix( /* args */ );
    }

  static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::D::double_create_CompRow_Matrix( /* args */ );
    }

  static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::D::double_create_Dense_Matrix( /* args */ );
    }
};


/* The specializations for Teuchos::as<> for New_Solver::complex and
 * New_Solver::doublecomplex are provided in Amesos2_NewSolver_Type.hpp
 */
template <>
struct FunctionMap<NewSolver,std::complex<float> >
{
  static void solve( /* args */ )
    {
      New_Solver::C::complex_solve( /* args */ );
    }

  static void symbolic_fact( /* args */ )
    {
      New_Solver::C::complex_symbolic_fact( /* args */ );
    }

  static void numeric_fact( /* args */ )
    {
      New_Solver::C::complex_numeric_fact( /* args */ );
    }

  static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_CompCol_Matrix( /* args */ );
    }

  static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_CompRow_Matrix( /* args */ );
    }

  static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::C::complex_create_Dense_Matrix( /* args */ );
    }
};


template <>
struct FunctionMap<NewSolver,std::complex<double> >
{
  static void solve( /* args */ )
    {
      New_Solver::Z::doublecomplex_solve( /* args */ );
    }

  static void symbolic_fact( /* args */ )
    {
      New_Solver::Z::doublecomplex_symbolic_fact( /* args */ );
    }

  static void numeric_fact( /* args */ )
    {
      New_Solver::Z::doublecomplex_numeric_fact( /* args */ );
    }

  static void create_CompCol_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_CompCol_Matrix( /* args */ );
    }

  static void create_CompRow_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_CompRow_Matrix( /* args */ );
    }

  static void create_Dense_Matrix( /* args */ )
    {
      New_Solver::Z::doublecomplex_create_Dense_Matrix( /* args */ );
    }
};

/*
 * etc. for as many types as need supporting
 */


} // end namespace Amesos

#endif  // AMESOS2_NEWSOLVER_FUNCTIONMAP_HPP
