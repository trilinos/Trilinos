// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_UNIT_TEST_HELPERS_HPP
#define STOKHOS_UNIT_TEST_HELPERS_HPP

#include "Stokhos_OrthogPolyApprox.hpp"
#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_FancyOStream.hpp"

namespace Stokhos {

  template<class OrdinalType, class ValueType>
  bool comparePCEs(const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>& a1,
                   const std::string& a1_name,
                   const Stokhos::OrthogPolyApprox<OrdinalType,ValueType>&a2,
                   const std::string& a2_name,
                   const ValueType& rel_tol, const ValueType& abs_tol,
                   Teuchos::FancyOStream& out)
  {
    bool success = true;

    out << "Comparing " << a1_name << " == " << a2_name << " ... ";

    const OrdinalType n = a1.size();

    // Compare sizes
    if (a2.size() != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == "
          << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }

    // Compare elements
    for( OrdinalType i = 0; i < n; ++i ) {
      ValueType nrm = std::sqrt(a1.basis()->norm_squared(i));
      ValueType err = std::abs(a1[i] - a2[i]) / nrm;
      ValueType tol =
        abs_tol + rel_tol*std::max(std::abs(a1[i]),std::abs(a2[i]))/nrm;
      if (err  > tol) {
        out
          <<"\nError, relErr("<<a1_name<<"["<<i<<"],"
          <<a2_name<<"["<<i<<"]) = relErr("<<a1[i]<<","<<a2[i]<<") = "
          <<err<<" <= tol = "<<tol<<": failed!\n";
        success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }
    else {
      out << std::endl
          << a1_name << " = " << a1 << std::endl
          << a2_name << " = " << a2 << std::endl;
    }

    return success;
  }

  template<class ValueType>
  bool compareValues(const ValueType& a1,
                     const std::string& a1_name,
                     const ValueType&a2,
                     const std::string& a2_name,
                     const ValueType& rel_tol, const ValueType& abs_tol,
                     Teuchos::FancyOStream& out)
  {
    bool success = true;

    ValueType err = std::abs(a1 - a2);
    ValueType tol = abs_tol + rel_tol*std::max(std::abs(a1),std::abs(a2));
    if (err  > tol) {
      out << "\nError, relErr(" << a1_name <<","
          << a2_name << ") = relErr(" << a1 <<"," << a2 <<") = "
          << err << " <= tol = " << tol << ": failed!\n";
      success = false;
    }

    return success;
  }

  template<class ValueType, class VectorType1, class VectorType2>
  bool compareVecs(const VectorType1& a1,
                   const std::string& a1_name,
                   const VectorType2&a2,
                   const std::string& a2_name,
                   const ValueType& rel_tol, const ValueType& abs_tol,
                   Teuchos::FancyOStream& out)
  {
    using value_t_1 = typename VectorType1::value_type;
    using value_t_2 = typename VectorType2::value_type;
    static_assert(std::is_same<value_t_1,value_t_2>::value,"Inconsistent types.");

    bool success = true;

    out << "Comparing " << a1_name << " == " << a2_name << " ... ";

    const int n = a1.size();

    // Compare sizes
    if (a2.size() != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == "
          << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }

    // Compare elements
    for( int i = 0; i < n; ++i ) {
      ValueType err = abs(a1.coeff(i) - a2.coeff(i));
      ValueType tol =
        abs_tol + rel_tol*std::max(abs(a1.fastAccessCoeff(i)),
                                   abs(a2.fastAccessCoeff(i)));
      if ( (err  > tol) && (err != 0.) ) {
        out
          <<"\nError, relErr("<<a1_name<<"["<<i<<"],"
          <<a2_name<<"["<<i<<"]) = relErr("<<a1.coeff(i)<<","<<a2.coeff(i)
          <<") = "<<err<<" <= tol = "<<tol<<": failed!\n";
        success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }
    else {
      out << std::endl
          << a1_name << " = " << a1 << std::endl
          << a2_name << " = " << a2 << std::endl;
    }

    return success;
  }

  template<class Array1, class Array2, class ValueType>
  bool compareArrays(const Array1& a1, const std::string& a1_name,
                     const Array2& a2, const std::string& a2_name,
                     const ValueType& rel_tol,
                     const ValueType& abs_tol,
                     Teuchos::FancyOStream& out)
  {
    using Teuchos::as;
    bool success = true;

    out << "Comparing " << a1_name << " == " << a2_name << " ... ";

    const int n = a1.size();

    // Compare sizes
    if (as<int>(a2.size()) != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == "
          << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }

    // Compare elements
    for( int i = 0; i < n; ++i ) {
      ValueType err = std::abs(a1[i] - a2[i]);
      ValueType tol =
        abs_tol + rel_tol*std::max(std::abs(a1[i]),std::abs(a2[i]));
      if (err > tol) {
        out << "\nError, relErr(" << a1_name << "[" << i << "]," << a2_name
            << "[" << i << "]) = relErr(" << a1[i] << "," <<a2[i] <<") = "
            << err << " <= tol = " << tol << ": failed!\n";
        success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }

    return success;
  }

  template<class ordinal_type, class scalar_type>
  bool compareSDM(
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& a1,
    const std::string& a1_name,
    const Teuchos::SerialDenseMatrix<ordinal_type,scalar_type>& a2,
    const std::string& a2_name,
    const scalar_type& rel_tol,
    const scalar_type& abs_tol,
    Teuchos::FancyOStream& out)
  {
    using Teuchos::as;
    bool success = true;

    out << "Comparing " << a1_name << " == " << a2_name << " ... ";

    const ordinal_type m = a1.numRows();
    const ordinal_type n = a1.numCols();

    // Compare number of rows
    if (a2.numRows() != m) {
      out << "\nError, "<<a1_name<<".numRows() = "<<a1.numRows()<<" == "
          << a2_name<<".numRows() = "<<a2.numRows()<<" : failed!\n";
      return false;
    }

    // Compare number of columnns
    if (a2.numCols() != n) {
      out << "\nError, "<<a1_name<<".numCols() = "<<a1.numCols()<<" == "
          << a2_name<<".numCols() = "<<a2.numCols()<<" : failed!\n";
      return false;
    }

    // Compare elements
    for (ordinal_type i=0; i<m; i++) {
      for (ordinal_type j=0; j<n; j++) {
        scalar_type err = std::abs(a1(i,j) - a2(i,j));
        scalar_type tol =
          abs_tol + rel_tol*std::max(std::abs(a1(i,j)),std::abs(a2(i,j)));
        if (err > tol) {
          out << "\nError, relErr(" << a1_name << "(" << i << "," << j << "), "
              << a2_name << "(" << i << "," << j << ")) = relErr("
              << a1(i,j) << ", " <<a2(i,j) <<") = "
              << err << " <= tol = " << tol << ": failed!\n";
          success = false;
        }
      }
    }
    if (success) {
      out << "passed\n";
    }

    return success;
  }

  template<class ordinal_type, class scalar_type>
  bool compareSparse3Tensor(
    const Stokhos::Sparse3Tensor<ordinal_type,scalar_type>& Cijk1,
    const std::string& cijk1_name,
    const Stokhos::Sparse3Tensor<ordinal_type,scalar_type>& Cijk2,
    const std::string& cijk2_name,
    const scalar_type& rel_tol,
    const scalar_type& abs_tol,
    Teuchos::FancyOStream& out)
  {
    typedef Stokhos::Sparse3Tensor<ordinal_type, scalar_type> Cijk_type;
    bool success = true;

    out << "Comparing " << cijk1_name << " == " << cijk2_name << " ... ";

     // Check number of nonzeros
    TEUCHOS_TEST_EQUALITY(Cijk1.num_entries(), Cijk2.num_entries(), out,
                          success);

    // Check entries
    for (typename Cijk_type::k_iterator k_it=Cijk2.k_begin();
         k_it!=Cijk2.k_end(); ++k_it) {
      ordinal_type k = Stokhos::index(k_it);
      for (typename Cijk_type::kj_iterator j_it = Cijk2.j_begin(k_it);
           j_it != Cijk2.j_end(k_it); ++j_it) {
        ordinal_type j = Stokhos::index(j_it);
        for (typename Cijk_type::kji_iterator i_it = Cijk2.i_begin(j_it);
             i_it != Cijk2.i_end(j_it); ++i_it) {
          ordinal_type i = Stokhos::index(i_it);
          scalar_type c1 = Cijk1.getValue(i,j,k);
          scalar_type c2 = Stokhos::value(i_it);

          double tol = abs_tol + c2*rel_tol;
          double err = std::abs(c1-c2);
          bool s = err < tol;
          if (!s) {
            out << std::endl
                << "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
                << " = " << "rel_err( " << c1 << ", " << c2 << " ) = " << err
                << " <= " << tol << " : ";
            if (s) out << "Passed.";
            else
              out << "Failed!";
            out << std::endl;
          }
          success = success && s;
        }
      }
    }

    return success;
  }

  template<class ordinal_type, class scalar_type>
  bool testSparse3Tensor(
    const Stokhos::Sparse3Tensor<ordinal_type,scalar_type>& Cijk,
    const Stokhos::ProductBasis<ordinal_type,scalar_type>& basis,
    const scalar_type& sparse_tol,
    const scalar_type& rel_tol,
    const scalar_type& abs_tol,
    Teuchos::FancyOStream& out,
    bool linear = false)
  {
    ordinal_type d = basis.dimension();
    ordinal_type sz = basis.size();
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,scalar_type> > > bases = basis.getCoordinateBases();
    Stokhos::Sparse3Tensor<ordinal_type,scalar_type> Cijk2;
    Teuchos::Array< Teuchos::RCP<Stokhos::Dense3Tensor<ordinal_type,scalar_type> > > Cijk_1d(d);
    for (ordinal_type i=0; i<d; i++)
      Cijk_1d[i] = bases[i]->computeTripleProductTensor();
    for (ordinal_type j=0; j<sz; j++) {
      Stokhos::MultiIndex<ordinal_type> terms_j = basis.term(j);
      for (ordinal_type i=0; i<sz; i++) {
        Stokhos::MultiIndex<ordinal_type> terms_i = basis.term(i);
        for (ordinal_type k=0; k<sz; k++) {
          Stokhos::MultiIndex<ordinal_type> terms_k = basis.term(k);
          if (!linear || terms_k.order() <= 1) {
            scalar_type c = 1.0;
            for (ordinal_type l=0; l<d; l++)
              c *= (*Cijk_1d[l])(terms_i[l],terms_j[l],terms_k[l]);
            if (std::abs(c/basis.norm_squared(i)) > sparse_tol)
              Cijk2.add_term(i,j,k,c);
          }
        }
      }
    }
    Cijk2.fillComplete();

    // std::cout << std::endl << Cijk << std::endl;
    // std::cout << std::endl << Cijk2 << std::endl;

    bool success = compareSparse3Tensor(Cijk, "Cijk", Cijk2, "Cijk2",
                                        rel_tol, abs_tol, out);

    return success;
  }

  template <typename operator_type1, typename operator_type2,
            typename scalar_type>
  bool testPseudoSpectralPoints(const operator_type1& op1,
                                const operator_type2& op2,
                                const scalar_type& rel_tol,
                                const scalar_type& abs_tol,
                                Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type1::const_set_iterator point_iterator_type;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(op1.point_size(), op2.point_size(), out, success);
    if (!success) return false;

    // Compare points
    point_iterator_type pi1 = op1.set_begin();
    point_iterator_type pi2 = op2.set_begin();
    while (pi1 != op1.set_end() || pi2 != op2.set_end()) {
      std::stringstream ss1, ss2;
      ss1 << pi1->first;
      ss2 << pi2->first;
      success = success &&
        Stokhos::compareArrays(pi1->first, ss1.str(),
                               pi2->first, ss2.str(),
                               rel_tol, abs_tol, out);

      std::stringstream ss3, ss4;
      ss3 << pi1->second.first;
      ss4 << pi2->second.first;
      success = success &&
        Stokhos::compareValues(pi1->second.first, ss3.str(),
                               pi2->second.first, ss4.str(),
                               rel_tol, abs_tol, out);
      ++pi1;
      ++pi2;
    }

    return success;
  }

  template <typename operator_type1, typename operator_type2,
            typename func_type1, typename func_type2,
            typename scalar_type>
  bool testPseudoSpectralApply(const operator_type1& op1,
                               const operator_type2& op2,
                               const func_type1& func1,
                               const func_type2& func2,
                               const scalar_type& rel_tol,
                               const scalar_type& abs_tol,
                               Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type1::ordinal_type ordinal_type;
    typedef typename operator_type1::value_type value_type;
    typedef typename operator_type1::const_iterator point_iterator_type;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(op1.coeff_size(), op2.coeff_size(), out, success);
    if (!success)
      return false;

    ordinal_type point_sz1 = op1.point_size();
    ordinal_type point_sz2 = op2.point_size();
    ordinal_type coeff_sz = op1.coeff_size();

    // Evaluate function at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(point_sz1,2);
    ordinal_type idx = 0;
    for (point_iterator_type pi = op1.begin(); pi != op1.end(); ++pi) {
      f(idx,0) = func1(*pi);
      f(idx,1) = func2(*pi);
      ++idx;
    }

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f2(point_sz2,2);
    idx = 0;
    for (point_iterator_type pi = op2.begin(); pi != op2.end(); ++pi) {
      f2(idx,0) = func1(*pi);
      f2(idx,1) = func2(*pi);
      ++idx;
    }

    // Compare function evaluations
    if (point_sz1 == point_sz2)
      success = success &&
        Stokhos::compareSDM(f, "f", f2, "f2", rel_tol, abs_tol, out);

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(coeff_sz,2);
    op1.transformQP2PCE(1.0, f, x, 0.0);

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x2(coeff_sz,2);
    op2.transformQP2PCE(1.0, f2, x2, 0.0);

    // Compare PCE coefficients
    success = success &&
      Stokhos::compareSDM(x, "x", x2, "x2", rel_tol, abs_tol, out);

    return success;
  }

  template <typename operator_type1, typename operator_type2,
            typename func_type1, typename func_type2,
            typename scalar_type>
  bool testPseudoSpectralApplyTrans(const operator_type1& op1,
                                    const operator_type2& op2,
                                    const func_type1& func1,
                                    const func_type2& func2,
                                    const scalar_type& rel_tol,
                                    const scalar_type& abs_tol,
                                    Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type1::ordinal_type ordinal_type;
    typedef typename operator_type1::value_type value_type;
    typedef typename operator_type1::const_iterator point_iterator_type;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(op1.coeff_size(), op2.coeff_size(), out, success);
    if (!success)
      return false;

    ordinal_type point_sz1 = op1.point_size();
    ordinal_type point_sz2 = op2.point_size();
    ordinal_type coeff_sz = op1.coeff_size();

    // Evaluate function at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(2,point_sz1);
    ordinal_type idx = 0;
    for (point_iterator_type pi = op1.begin(); pi != op1.end(); ++pi) {
      f(0,idx) = func1(*pi);
      f(1,idx) = func2(*pi);
      ++idx;
    }

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f2(2,point_sz2);
    idx = 0;
    for (point_iterator_type pi = op2.begin(); pi != op2.end(); ++pi) {
      f2(0,idx) = func1(*pi);
      f2(1,idx) = func2(*pi);
      ++idx;
    }

    // Compare function evaluations
    if (point_sz1 == point_sz2)
      success = success &&
        Stokhos::compareSDM(f, "f", f2, "f2", rel_tol, abs_tol, out);

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(2,coeff_sz);
    op1.transformQP2PCE(1.0, f, x, 0.0, true);

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x2(2,coeff_sz);
    op2.transformQP2PCE(1.0, f2, x2, 0.0, true);

    // Compare PCE coefficients
    success = success &&
      Stokhos::compareSDM(x, "x", x2, "x2", rel_tol, abs_tol, out);

    return success;
  }

  template <typename basis_type, typename operator_type, typename scalar_type>
  bool testPseudoSpectralDiscreteOrthogonality(const basis_type& basis,
                                               const operator_type& op,
                                               const scalar_type& rel_tol,
                                               const scalar_type& abs_tol,
                                               Teuchos::FancyOStream& out) {
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    ordinal_type coeff_sz = op.coeff_size();
    ordinal_type point_sz = op.point_size();
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> eye(coeff_sz,coeff_sz);
    eye.putScalar(0.0);
    for (ordinal_type i=0; i<coeff_sz; ++i)
      eye(i,i) = 1.0;

    // Map PCE coefficients to quad values
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(point_sz,coeff_sz);
    op.transformPCE2QP(1.0, eye, f, 0.0);

    // Map quad values to PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(coeff_sz,coeff_sz);
    op.transformQP2PCE(1.0, f, x, 0.0);

    // Subtract identity
    for (ordinal_type i=0; i<coeff_sz; ++i)
      x(i,i) -= 1.0;

    // Expected answer, which is zero
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> z(coeff_sz,coeff_sz);
    z.putScalar(0.0);

    out << "Discrete orthogonality error = " << x.normInf() << std::endl;

    // Compare PCE coefficients
    bool success = Stokhos::compareSDM(x, "x", z, "zero", rel_tol, abs_tol, out);

    return success;
  }

  /*
  template <typename basis_type, typename operator_type, typename scalar_type>
  bool testPseudoSpectralDiscreteOrthogonality(const basis_type& basis,
                                               const operator_type& op,
                                               const scalar_type& rel_tol,
                                               const scalar_type& abs_tol,
                                               Teuchos::FancyOStream& out) {
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    // Evaluate full basis at all quadrature points
    ordinal_type coeff_sz = op.coeff_size();
    ordinal_type point_sz = op.point_size();
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(point_sz,coeff_sz);
    Teuchos::Array<value_type> vals(coeff_sz);
    ordinal_type i = 0;
    typename operator_type::const_iterator pi = op.begin();
    for (; pi != op.end(); ++pi) {
      basis.evaluateBases(*pi, vals);
      for (ordinal_type j=0; j<coeff_sz; ++j)
        f(i,j) = vals[j];

      ++i;
    }

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(coeff_sz,coeff_sz);
    op.transformQP2PCE(1.0, f, x, 0.0);

    // Subtract identity
    for (ordinal_type i=0; i<coeff_sz; ++i)
      x(i,i) -= 1.0;

    // Expected answer, which is zero
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> z(coeff_sz,coeff_sz);
    z.putScalar(0.0);

    out << "Discrete orthogonality error = " << x.normInf() << std::endl;

    // Compare PCE coefficients
    bool success = Stokhos::compareSDM(x, "x", z, "zero", 1e-14, 1e-14, out);

    return success;
  }
  */
}

#endif // STOKHOS_UNIT_TEST_HELPERS_HPP
