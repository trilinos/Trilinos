#ifndef ROL_teuchosarray
#define ROL_teuchosarray

#include "Teuchos_Array.hpp"
#include <iostream>

// Defines the vector space used for optimization.
namespace ROL{
    struct teuchosarrayVS {
        typedef double Real;
        typedef Teuchos::Array <double> Vector;

        // y <- x (Shallow.  No memory allocation.)
        static void copy(const Vector& x, Vector& y){
            for(unsigned int i=0;i<x.size();i++){
                y[i]=x[i];
            }
        }

        // x <- alpha * x
        static void scal(const Real& alpha, Vector& x){
            for(unsigned int i=0;i<x.size();i++){
                x[i]=alpha*x[i];
            }
        }

        // y <- alpha * x + y
        static void axpy(const Real& alpha, const Vector& x, Vector& y){
            for(unsigned int i=0;i<x.size();i++){
                y[i]=alpha*x[i]+y[i];
            }
        }

        // innr <- <x,y>
        static Real innr(const Vector& x,const Vector& y){
            Real z=0;
            for(unsigned int i=0;i<x.size();i++){
                z+=x[i]*y[i];
            }
            return z;
        }

        // Memory allocation and size setting
        static void init(const Vector& x, Vector& y){
            y.resize(x.size());
        }

        // Prints out diagnostic information
        static void print(std::string msg){
            std::cout << msg;
        }

        // Prints out error information
        static void error(std::string msg){
            std::cout << msg;
        }
    };
}

#endif
