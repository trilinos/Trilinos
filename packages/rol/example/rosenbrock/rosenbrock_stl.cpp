#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include "ROL_Algorithms.hpp"
#include "ROL_stl.hpp"

typedef ROL::stlVS MyVS;

double sq(double x){
    return x*x; 
}

// Objective function for Rosenbrock
class RosenObjective : public ROL::Functional <MyVS> {
public:
    double operator () (const std::vector<double>& x) const{
        return sq(1.-x[0])+100.*sq(x[1]-sq(x[0]));
    }
};

// Gradient operator for Rosenbrock
class RosenGradient : public ROL::Operator <MyVS,MyVS> {
public:
    void operator () (
        const std::vector<double>& x,
        std::vector<double>& g
    ) const{
        g[0]=-400*x[0]*(x[1]-sq(x[0]))-2*(1-x[0]);
        g[1]=200*(x[1]-sq(x[0]));
    }
};

// Hessian operator for Rosenbrock
class RosenHessian: public ROL::Operator <MyVS,MyVS> {
private:
    // Current optimization iterate
    const std::vector <double>& x;
public:
    RosenHessian(const std::vector <double>& x_) : x(x_) {};

    void operator () (
        const std::vector<double>& eta,
        std::vector<double>& Heta 
    ) const{
        Heta[0]= (1200*sq(x[0])-400*x[1]+2)*eta[0]
            -400*x[0]*eta[1];

        Heta[1]= -400*x[0]*eta[0] + 200*eta[1];
    }
};

int main(){

    // Generate an initial guess for Rosenbrock
    std::vector <double> x(2);
    x[0]=-1.2; x[1]=1.;

    // Create a direction for the finite difference tests
    std::vector <double> eta(2);
    eta[0]=-.5; eta[1]=.5;
    
    // Create a state and setup the problem
    ROL::core<MyVS>::State state(x);
    state.H_type=ROL::External_t;
    state.algorithm_class=ROL::TrustRegion;
    state.eps_g=1e-10;
    state.eps_s=1e-10;
    state.iter_max=200;
    state.eps_krylov=1e-8;

    // Create a function, gradient, and Hessian for this problem
    RosenObjective F;
    RosenGradient G;
    // Make sure to link the Hessian to the current optimization
    // iterate
    RosenHessian H(*(state.u.begin()));
   
    // Do a finite difference test for the gradient and Hessian
    std::cout << "Gradient finite difference check" << std::endl;
    ROL::derivativeCheck<MyVS>(F,G,x,eta);

    std::cout << "Hessian finite difference check" << std::endl;
    ROL::derivativeCheck<MyVS,MyVS>(G,H,x,eta,x);

    // Optimize the problem
    ROL::core<MyVS>::getMin(state,F,G,H);

    // Print out the final answer
    const std::vector <double>& opt_x=*(state.u.begin());
    std::cout << "The optimal point is: (" << opt_x[0] << ','
        << opt_x[1] << ')' << std::endl;

    std::cout << "End Result: TEST PASSED\n";
}
