// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test ROL::LA::Matrix and ROL::LA::Vector compatibility with LAPACK
*/

#include <iostream>
#include <vector>
#include <cmath>

#include "ROL_Stream.hpp"
#include "ROL_LinearAlgebra.hpp"
#include "ROL_LAPACK.hpp"

int main(int argc, char* argv[]) {
    
    // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
    auto outStream = ROL::makeStreamPtr(std::cout, argc > 1);
    
    int errorFlag = 0;
    
    try {
        
        *outStream << "Testing ROL::LA::Matrix LAPACK compatibility..." << std::endl;
        
        // Test 1: Basic matrix properties
        *outStream << "\n=== Test 1: Basic Matrix Properties ===" << std::endl;
        
        ROL::LA::Matrix<double> A(3, 3);
        A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;
        A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;
        A(2,0) = 7.0; A(2,1) = 8.0; A(2,2) = 10.0; // Make it slightly non-singular
        
        *outStream << "Matrix dimensions: " << A.numRows() << "x" << A.numCols() << std::endl;
        
        // Test 2: Memory layout and data access
        *outStream << "\n=== Test 2: Memory Layout ===" << std::endl;
        
        double* data = A.values();
        *outStream << "Raw data pointer validity: " << (data != nullptr ? "OK" : "FAIL") << std::endl;
        
        if (data == nullptr) {
            errorFlag++;
        }
        
        // Check if we can access matrix elements correctly
        *outStream << "Element access test:" << std::endl;
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < 3; j++) {
                *outStream << "A(" << i << "," << j << ") = " << A(i,j) << std::endl;
            }
        }
        
        // Test 3: LAPACK GEEV eigenvalue computation
        *outStream << "\n=== Test 3: LAPACK GEEV Call ===" << std::endl;
        
        ROL::LAPACK<int, double> lapack;
        ROL::LA::Matrix<double> mymat(A);  // Copy for GEEV (it will be destroyed)
        
        char jobvl = 'N';
        char jobvr = 'N';
        int n = 3;
        
        std::vector<double> real(n, 0);
        std::vector<double> imag(n, 0);
        
        double* vl = nullptr;
        double* vr = nullptr;
        int ldvl = 1;
        int ldvr = 1;
        int lwork = 4*n;
        std::vector<double> work(lwork, 0);
        int info = 0;
        
        *outStream << "Matrix leading dimension (stride): " << mymat.stride() << std::endl;
        *outStream << "Expected leading dimension: " << n << std::endl;
        
        if (mymat.stride() != n) {
            *outStream << "WARNING: Stride mismatch! This could cause LAPACK issues." << std::endl;
        }
        
        // Call GEEV
        lapack.GEEV(jobvl, jobvr, n, &mymat(0,0), n, &real[0], &imag[0], vl, ldvl, vr, ldvr, &work[0], lwork, &info);
        
        *outStream << "GEEV return code (info): " << info << std::endl;
        if (info != 0) {
            *outStream << "ERROR: GEEV failed with info = " << info << std::endl;
            errorFlag++;
        }
        
        *outStream << "Eigenvalues:" << std::endl;
        for(int i = 0; i < n; i++) {
            if (std::abs(imag[i]) < 1e-12) {
                *outStream << "  λ[" << i << "] = " << real[i] << std::endl;
            } else {
                *outStream << "  λ[" << i << "] = " << real[i] << " + " << imag[i] << "i" << std::endl;
            }
        }
        
        // Test 4: Comparison with known values
        *outStream << "\n=== Test 4: Eigenvalue Validation ===" << std::endl;
        
        // For this specific 3x3 matrix, we expect real eigenvalues
        // Let's just check that they're reasonable (trace = sum of eigenvalues)
        double trace = A(0,0) + A(1,1) + A(2,2);  // Should be 1 + 5 + 10 = 16
        double eigenvalue_sum = real[0] + real[1] + real[2];
        
        *outStream << "Matrix trace: " << trace << std::endl;
        *outStream << "Sum of eigenvalues: " << eigenvalue_sum << std::endl;
        *outStream << "Difference: " << std::abs(trace - eigenvalue_sum) << std::endl;
        
        if (std::abs(trace - eigenvalue_sum) > 1e-10) {
            *outStream << "ERROR: Trace and eigenvalue sum differ significantly!" << std::endl;
            errorFlag++;
        }
        
        // Test 5: Matrix stride/memory layout verification
        *outStream << "\n=== Test 5: Memory Layout Verification ===" << std::endl;
        
        ROL::LA::Matrix<double> B(2, 3);
        B(0,0) = 1; B(0,1) = 2; B(0,2) = 3;
        B(1,0) = 4; B(1,1) = 5; B(1,2) = 6;
        
        *outStream << "2x3 Matrix B:" << std::endl;
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < 3; j++) {
                *outStream << B(i,j) << " ";
            }
            *outStream << std::endl;
        }
        
        double* bdata = B.values();
        *outStream << "Raw memory layout: ";
        for(int k = 0; k < 6; k++) {
            *outStream << bdata[k] << " ";
        }
        *outStream << std::endl;
        
        *outStream << "Stride: " << B.stride() << std::endl;
        
        // For column-major storage, we expect: 1, 4, 2, 5, 3, 6
        // For row-major storage, we expect: 1, 2, 3, 4, 5, 6
        bool is_column_major = (bdata[0] == 1 && bdata[1] == 4 && bdata[2] == 2);
        bool is_row_major = (bdata[0] == 1 && bdata[1] == 2 && bdata[2] == 3);
        
        *outStream << "Matrix appears to be: ";
        if (is_column_major) {
            *outStream << "Column-major" << std::endl;
        } else if (is_row_major) {
            *outStream << "Row-major" << std::endl;
        } else {
            *outStream << "Unknown layout!" << std::endl;
            errorFlag++;
        }
        
    }
    catch (std::exception& e) {
        *outStream << "Exception caught: " << e.what() << std::endl;
        errorFlag++;
    }
    
    if (errorFlag == 0) {
      std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    } else {
      std::cout << "\nEnd Result: TEST FAILED (errors: " << errorFlag << ")" << std::endl;
    }
    
    return 0;
}
