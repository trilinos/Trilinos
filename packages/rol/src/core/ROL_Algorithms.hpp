/* Included, licensed code
peopt .1 beta
Author: Joseph Young (josyoun@sandia.gov)
Copyright 2012 Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains certain
rights in this software.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
    * this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
    * notice, this list of conditions and the following disclaimer in the
    * documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ROL_H
#define ROL_H

#include<list>
#include<map>
#include<limits>
#include<cmath>
#include<sstream>
#include<iomanip>
#include<memory>

namespace ROL{
    template <typename VS>
    class DataStructures{
        // No copy constructor required
    // We really want this to be private, but AIX throws an error
    public:
        // Prevent allocation of this class
        explicit DataStructures();
        // No assignment operator required

        template<typename T>
        class Alloc : public std::allocator <T> {
        public:
            typedef std::allocator <T> base;
            typedef typename base::size_type        size_type;
            typedef typename base::difference_type  difference_type;
            typedef typename base::pointer          pointer;
            typedef typename base::const_pointer    const_pointer;
            typedef typename base::reference        reference;
            typedef typename base::const_reference  const_reference;
            typedef typename base::value_type       value_type;

            template<typename U>
            struct rebind {
                typedef Alloc <U> other;
            };

            Alloc() {}

            // In this case, we can not make the constructor explicit
            template<typename U>
            Alloc(const Alloc<U>&) {}

            void construct(pointer p, const_reference t){
               T* q=new( (void*)p ) T();
               VS::init(t,*q);
               VS::copy(t,*q);
            }
        };

    public:
        typedef std::list<typename VS::Vector,Alloc<typename VS::Vector> > List;
        typedef std::map<typename VS::Vector,Alloc<typename VS::Vector> > Map;
    };

    // A simple operator specification 
    template <class Domain, class Codomain>
    class Operator {
    public:
        // Basic application
        virtual void operator () (
            const typename Domain::Vector& x,
            typename Codomain::Vector &y) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~Operator() {}
    };

    // A simple functional interface
    template <class Domain>
    class Functional {
    public:
        // Basic application
        virtual typename Domain::Real operator ()
            (const typename Domain::Vector& x) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~Functional() {}
    };

    // A combination of the operator and functional functionalities
    template <class Domain,class Codomain>
    class OperatorFunctional : public Operator <Domain,Codomain> {
    public:
        // Basic application
        virtual void operator () (
            const typename Domain::Vector& x,
            typename Codomain::Vector& y,
            typename Domain::Real& obj_val) const = 0;

        // Allow a derived class to deallocate memory 
        virtual ~OperatorFunctional() {}
    };
        
    // Operator possessing two derivatives and adjoints of those derivatives
    template <class Domain, class Codomain>
    class DiffOperator {
    public:
        // Basic application
        virtual std::auto_ptr < Operator<Domain,Codomain> >
            f(const int i) const=0;

        // Derivative of the operator
        virtual std::auto_ptr < Operator<Domain,Codomain> >
            fp(const int i,const typename Domain::Vector& x) const =0;

        // Derivative adjoint
        virtual std::auto_ptr < Operator<Codomain,Domain> >
            fps(const int i,const typename Domain::Vector& x) const =0;
        
        // Second derivative in the direction eta, adjoint, in the
        // direction xi
        virtual std::auto_ptr < Operator<Domain,Domain> >
            fpps(
                const int i,
                const typename Domain::Vector& x,
                const typename Codomain::Vector& xi) const =0;

        // Gets the maximum index for each of the above operators
        virtual int max_index() const=0;

        // Allow the collection of operators to free memory
        virtual ~DiffOperator() {}
    };

    // Performs a 4-point finite difference directional derivative on
    // a scalar valued function.
    template <class Domain>
    static typename Domain::Real directionalDerivative(
        const Functional<Domain>& F,
        const typename Domain::Vector& u,
        const typename Domain::Vector& eta,
        const typename Domain::Real& epsilon
    ){
            // Create some type short-cuts
        typedef typename Domain::Vector Vector;
        typedef typename Domain::Real Real;

        // Create a single work element
        Vector work; Domain::init(u,work);

        // F(u+eps s)
        Domain::copy(u,work);
        Domain::axpy(epsilon,eta,work);
        Real obj_upes=F(work);

        // F(u-eps s)
        Domain::copy(u,work);
        Domain::axpy(-epsilon,eta,work);
        Real obj_umes=F(work);

        // F(u+2 eps s)
        Domain::copy(u,work);
        Domain::axpy(Real(2.*epsilon),eta,work);
        Real obj_up2es=F(work);

        // F(u-2 eps s)
        Domain::copy(u,work);
        Domain::axpy(Real(-2.*epsilon),eta,work);
        Real obj_um2es=F(work);

        // Calculate the directional derivative and return it
        Real dd=(obj_um2es-Real(8.)*obj_umes+Real(8.)*obj_upes-obj_up2es)
            /(Real(12.)*epsilon);
        return dd;
    }

    // Performs a 4-point finite difference directional derivative on
    // a vector valued function.
    template <class Domain,class Codomain>
    static void directionalDerivative(
        const Operator<Domain,Codomain>& F,
        const typename Domain::Vector& u,
        const typename Domain::Vector& eta,
        const typename Domain::Real& epsilon,
        typename Codomain::Vector& dd
    ){
            // Create some type shortcuts
        typedef typename Domain::Real Domain_Real;
        typedef typename Codomain::Real Codomain_Real;

        // Zero out the directional derivative
        Codomain::scal(Codomain_Real(0.),dd);

        // Create two work elements 
        typename Domain::Vector work; Domain::init(u,work);
        typename Codomain::Vector work2; Codomain::init(dd,work2);

        // F(u+eps s)
        Domain::copy(u,work);
        Domain::axpy(epsilon,eta,work);
        F(work,work2);
        Codomain::axpy(Codomain_Real(8.),work2,dd);

        // F(u-eps s)
        Domain::copy(u,work);
        Domain::axpy(-epsilon,eta,work);
        F(work,work2);
        Codomain::axpy(Codomain_Real(-8.),work2,dd);

        // F(u+2 eps s)
        Domain::copy(u,work);
        Domain::axpy(Domain_Real(2.)*epsilon,eta,work);
        F(work,work2);
        Codomain::axpy(Codomain_Real(-1.),work2,dd);

        // F(u-2 eps s)
        Domain::copy(u,work);
        Domain::axpy(Domain_Real(-2.)*epsilon,eta,work);
        F(work,work2);
        Codomain::axpy(Codomain_Real(1.),work2,dd);

        // Finish the finite difference calculation 
        Codomain::scal(Codomain_Real(1.)/(Codomain_Real(12.)*epsilon),dd);
    }

    // Performs a finite difference test between two operators F and G where
    // F is scalar valued.  This helps verify that G is the derivative of F.
    // In other words, we check that grad(F)(u)=G(u).
    template <class Domain>
    static void derivativeCheck(
        Functional<Domain>& F,
        Operator<Domain,Domain>& G,
        const typename Domain::Vector& u,
        const typename Domain::Vector& eta
    ) {
        // Create some type shortcuts
        typedef typename Domain::Vector Vector;
        typedef typename Domain::Real Real;

        // Calculate the gradient at the point u
        Vector g; Domain::init(u,g);
        G(u,g);

        // Begin by calculating the directional derivative via the gradient
        Real dd_grad=Domain::innr(g,eta);

        // Compute an ensemble of finite difference tests in a linear manner
        Domain::print("Finite difference test of a gradient.\n");
        for(int i=-2;i<=5;i++){
            typename Domain::Real epsilon=pow(Real(.1),i);
            typename Domain::Real dd=
                directionalDerivative <Domain> (F,u,eta,epsilon);

            std::stringstream ss;
            if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
            else ss << "The relative difference (1e-" << i << "): ";
            ss << std::scientific << std::setprecision(16)
                << fabs(dd_grad-dd)/(Real(1e-16)+fabs(dd_grad))<< std::endl;
            Domain::print(ss.str());
        }
        Domain::print("\n");
    }

    // Performs a finite difference test between two operators F and G where
    // F is vector valued.  This helps verify that G is the derivative of F.
    // In other words, we check that F'(u)eta = G(u)eta.
    template <class Domain,class Codomain>
    static void derivativeCheck(
        Operator<Domain,Codomain>& F,
        Operator<Domain,Codomain>& G,
        const typename Domain::Vector& u,
        const typename Domain::Vector& eta,
        const typename Domain::Vector& y
    ) {
        // Create some type shortcuts
        typedef typename Codomain::Real Codomain_Real;
        typedef typename Domain::Real Domain_Real;

        // Create a single work element
        typename Codomain::Vector work; Codomain::init(y,work);

        // Calculate the operator G in the direction eta.  Typically, this
        // is something like a Hessian-vector product.  Here, we assume that
        // the base of G is set to u.
        typename Codomain::Vector G_eta; Codomain::init(y,G_eta);
        G(eta,G_eta);

        // Compute an ensemble of finite difference tests in a linear manner
        Domain::print("Finite difference test of an operator-vector "
            "product.\n");
        for(int i=-2;i<=5;i++){

            // Calculate the directional derivative
            Domain_Real epsilon=pow(Domain_Real(.1),i);
            directionalDerivative <Domain,Codomain> (F,u,eta,epsilon,work);

            // Determine the residual.  Store in work.
            Codomain::axpy(Codomain_Real(-1.),G_eta,work);

            // Determine the relative error
            Codomain_Real rel_err=sqrt(Codomain::innr(work,work))
                /(Codomain_Real(1e-16)+sqrt(Codomain::innr(G_eta,G_eta)));

            // Print out the differences
            std::stringstream ss;
            if(i<0) ss << "The relative difference (1e+" << -i <<  "): ";
            else ss << "The relative difference (1e-" << i << "): ";
            ss << std::scientific << std::setprecision(16)
                << rel_err << std::endl;
            Domain::print(ss.str());
        }
        Domain::print("\n");
    }

    // Which algorithm class do we use
    enum AlgorithmClass{
        TrustRegion,            // Trust-Region algorithms
        LineSearch              // Line-search algorithms
    };

    // Reasons why we stop the algorithm
    enum StoppingCondition{
        NotConverged,            // Algorithm did not converge
        RelativeGradientSmall,   // Relative gradient was sufficiently small
        RelativeStepSmall,       // Relative change in the step is small
        MaxItersExceeded         // Maximum number of iterations exceeded
    };

    // Reasons we stop the Krylov method
    enum KrylovStop{
        NegativeCurvature,        // Negative curvature detected
        RelativeErrorSmall,       // Relative error is small
        MaxKrylovItersExceeded,   // Maximum number of iterations exceeded
        TrustRegionViolated       // Trust-region radius violated
    };

    // Various operators for both Hessian approximations and preconditioners
    enum Operators{
        Identity_t,        // Identity approximation
        ScaledIdentity_t,  // Identity approximation
        BFGS_t,            // BFGS approximation
        InvBFGS_t,         // Inverse BFGS approximation
        SR1_t,             // SR1 approximation
        InvSR1_t,          // Inverse SR1 approximation
        External_t           // An external operator provided by the user
    };

    // Different kinds of search directions 
    enum LineSearchDirection{
        SteepestDescent_t,        // SteepestDescent 
        FletcherReeves_t,         // Fletcher-Reeves CG
        PolakRibiere_t,           // Polak-Ribiere CG
        HestenesStiefel_t,        // HestenesStiefel CG
        LimitedMemoryBFGS_t,      // Limited-memory BFGS 
        NewtonCG_t                // Newton-CG
    };

    enum LineSearchKind{
        Brents_t,         // Brent's minimization
        GoldenSection_t,  // Golden-section search 
        BackTracking_t,   // BackTracking search 
        TwoPointA_t,      // Barzilai and Borwein's method A
        TwoPointB_t       // Barzilai and Borwein's method B
    };

    // The core routines for rol 
    template <typename VS>
    struct core{
        // Create some type shortcuts 
        typedef typename VS::Real Real;
        typedef typename VS::Vector Vector;
        typedef typename DataStructures <VS>::List List;
        typedef typename DataStructures <VS>::List::iterator ListIterator;
        typedef typename DataStructures <VS>::List::const_iterator
            ListConstIterator;

        // Pieces of the state required for all optimization routines        
        struct State{
        
            // ------------- GENERIC ------------- 

            // Tolerance for the gradient stopping condition
            Real eps_g;

            // Tolerance for the step length stopping criteria
            Real eps_s;

            // Number of control objects to store in a quasi-Newton method
            unsigned int stored_history;

            // Number of failed iterations before we reset the history for
            // quasi-Newton methods
            int history_reset;

            // Current iteration
            int iter;

            // Maximum number of optimization iterations
            int iter_max;

            // Why we've stopped the optimization
            StoppingCondition opt_stop;

            // Current number of Krylov iterations taken
            int krylov_iter;

            // Maximum number of iterations in the Krylov method
            int krylov_iter_max;

            // Total number of Krylov iterations taken
            int krylov_iter_total;

            // Why the Krylov method was last stopped
            KrylovStop krylov_stop;

            // Relative error in the Krylov method
            Real krylov_rel_err;

            // Stopping tolerance for the Krylov method
            Real eps_krylov;

            // Algorithm class
            AlgorithmClass algorithm_class;

            // Preconditioner
            Operators Minv_type;

            // Hessian approximation
            Operators H_type;

            // Norm of the gradient
            Real norm_g;

            // Norm of a typical tradient
            Real norm_gtyp;

            // Norm of the trial step
            Real norm_s;

            // Norm of a typical trial step
            Real norm_styp;

            // Optimization variable 
            List u; 
            
            // Gradient 
            List g;
            
            // Trial step 
            List s;
            
            // Old optimization variable 
            List u_old; 
            
            // Old gradient 
            List g_old;
            
            // Old trial step 
            List s_old;

            // Contains the prior iteration information for the
            // quasi-Newton operators
            List oldY;
            List oldS;

            // Current objective value
            Real obj_u;

            // Objective value at the trial step
            Real obj_ups;

            // Amount of verbosity
            int verbose;
            
            // ------------- TRUST-REGION ------------- 

            // Trust region radius
            Real delta;

            // Maximum trust region radius
            Real delta_max;

            // Trust-region parameter for checking whether a step has been
            // accepted
            Real eta1;

            // Trust-region parameter for checking whether a step has been
            // accepted
            Real eta2;

            // Ratio between the predicted and actual reduction
            Real rho;

            // Number of rejected trust-region steps
            int rejected_trustregion;

            // ------------- LINE-SEARCH ------------- 

            // Line-search step length
            Real alpha;

            // Current number of iterations used in the line-search
            int linesearch_iter;

            // Maximum number of iterations used in the line-search
            int linesearch_iter_max;

            // Total number of line-search iterations computed
            int linesearch_iter_total;

            // Stopping tolerance for the line-search
            double eps_ls;

            // Search direction type
            LineSearchDirection dir;

            // Type of line-search 
            LineSearchKind kind;

            // Initialize the state without setting up any variables. 
            State(){ init(); };

            // Initialize the state with the initial guess of the variables x
            State(const Vector& x) {
                init();
                u.push_back(x); 
                g.push_back(x);
                s.push_back(x);
                u_old.push_back(x);
                g_old.push_back(x);
                s_old.push_back(x);
            }
        private:
            void init(){
                eps_g=Real(1e-6);
                eps_s=Real(1e-6);
                stored_history=0;
                history_reset=5;
                iter=1;
                iter_max=10;
                opt_stop=NotConverged;
                krylov_iter=0;
                krylov_iter_max=10;
                krylov_iter_total=0;
                krylov_stop=RelativeErrorSmall;
                krylov_rel_err=Real(std::numeric_limits<Real>::quiet_NaN());
                eps_krylov=Real(1e-2);
                algorithm_class=TrustRegion;
                Minv_type=Identity_t;
                H_type=Identity_t;
                norm_g=Real(std::numeric_limits<Real>::quiet_NaN());
                norm_gtyp=Real(std::numeric_limits<Real>::quiet_NaN());
                norm_s=Real(std::numeric_limits<Real>::quiet_NaN());
                norm_styp=Real(std::numeric_limits<Real>::quiet_NaN());
                obj_u=Real(std::numeric_limits<Real>::quiet_NaN());
                obj_ups=Real(std::numeric_limits<Real>::quiet_NaN());
                verbose=1;
                delta=Real(100.);
                delta_max=Real(100.);
                eta1=Real(.1);
                eta2=Real(.9);
                rho=Real(std::numeric_limits<Real>::quiet_NaN());
                rejected_trustregion=0;
                alpha=1.;
                linesearch_iter=0;
                linesearch_iter_max=5;
                linesearch_iter_total=0;
                eps_ls=Real(1e-2);
                dir=SteepestDescent_t;
                kind=GoldenSection_t;
            }
        };
    
        // A function that has free reign to manipulate or analyze the state.
        // This should be used cautiously.
        class StateManipulator {
            // No default assignment operator or copy constructor.  The user
            // can define one, but the optimization does not require it.
        public:
            // Application
            virtual void operator () (State& state) {};

            // Allow the derived class to deallocate memory
            virtual ~StateManipulator() {}
        };

        // The identity operator 
        class Identity : public Operator <VS,VS> {
        public:
            void operator () (const Vector& p, Vector& result) const{
                VS::copy(p,result);
            }
        };

        // The scaled identity Hessian approximation.  Specifically, use use
        // norm(g) / delta_max I.
        class ScaledIdentity : public Operator <VS,VS> {
        private:
            const Real& norm_g;
            const Real& delta_max;
        public:
            explicit ScaledIdentity(State& state)
                : norm_g(state.norm_g), delta_max(state.delta_max) {};

            void operator () (const Vector& p, Vector& result) const{
                VS::copy(p,result);
                VS::scal(norm_g/delta_max,result);
            }
        };

        // The BFGS Hessian approximation.  
        /* Note, the formula we normally see for BFGS denotes the inverse
            Hessian approximation.  This is not the inverse, but the true
            Hessian approximation. */ 
        class BFGS : public Operator <VS,VS> {
        private:
            const List& oldY;
            const List& oldS;
        public:
            explicit BFGS(const State& state)
                : oldY(state.oldY), oldS(state.oldS) {};

            // Operator interface
            /* It's not entirely clear to me what the best implementation for
                this method really is.  In the following implementation, we
                require an additional k work elements where k is the number of
                stored gradient and position differences.  It's possible to
                reduce this to 1 or 2, but we need to compute redundant
                information.  It's also possible to implementation the compact
                representation, see "Representations of quasi-Newton matrices
                and their use in limited memory methods" from Byrd, Nocedal,
                and Schnabel.  The problem with that algorithm is that is
                requires machinery such as linear system solves that we don't
                current have.  It also works much better with matrices or
                multivectors of data and we don't require the user to provide
                these abstractions. */
            void operator () (const Vector& p, Vector& result) const{

                // Check that the number of stored gradient and trial step
                // differences is the same.
                if(oldY.size() != oldS.size())
                    VS::error("In the BFGS Hessian approximation, the number "
                        "of stored gradient differences must equal the number "
                        "of stored trial step differences.");

                // Allocate memory for work
                List work(oldY.size(),p);

                // If we have no vectors in our history, we return the direction
                VS::copy(p,result);
                if(oldY.size() == 0) return;

                // As a safety check, insure that the inner product between all
                // the (s,y) pairs is positive
                ListConstIterator y0=oldY.begin();
                ListConstIterator s0=oldS.begin();
                while(y0!=oldY.end()){
                    Real inner_y_s=VS::innr(*y0++,*s0++);
                    if(inner_y_s<0)
                        VS::error("Detected a (s,y) pair in BFGS that possesed "
                        "a nonpositive inner product");
                }

                // Othwerwise, we copy all of the trial step differences into
                // the work space
                ListIterator Bisj_iter=work.begin();
                ListConstIterator sk_iter=oldS.begin();
                while(Bisj_iter!=work.end())
                    VS::copy((*sk_iter++),(*Bisj_iter++));

                // Keep track of the element Bisi
                ListConstIterator Bisi_iter=work.end(); Bisi_iter--;

                // Keep iterating until Bisi equals the first element in the
                // work list. This means we have computed B1s1, B2s2, ..., Bksk.
                Bisj_iter=work.begin();
                ListConstIterator si_iter=oldS.end(); si_iter--;
                ListConstIterator yi_iter=oldY.end(); yi_iter--;
                ListConstIterator sj_iter=oldS.begin();
                while(1){

                    // Create some reference to our iterators that are easier to
                    // work with
                    const Vector& si=*si_iter;
                    const Vector& yi=*yi_iter;
                    const Vector& Bisi=*Bisi_iter;

                    // Determine <Bisi,si>
                    Real inner_Bisi_si=VS::innr(Bisi,si);

                    // Determine <yi,si>
                    Real inner_yi_si=VS::innr(yi,si);

                    // Determine <si,Bip>
                    Real inner_si_Bip=VS::innr(si,result);

                    // Determine <yi,p>
                    Real inner_yi_p=VS::innr(yi,p);

                    // Determine -<si,Bip>/<Bisi,si> Bisi + Bip.  Store in Bip.
                    // This will become B_{i+1}p.
                    VS::axpy(-inner_si_Bip/inner_Bisi_si,Bisi,result);

                    // Determine <yi,p>/<yi,si> yi + w where we calculated w
                    // in the line above.  This completes the calculation of
                    // B_{i+1}p
                    VS::axpy(inner_yi_p/inner_yi_si,yi,result);

                    // Check whether or not we've calculated B_{i+1}p for the
                    // last time
                    if(Bisi_iter==work.begin()) break;

                    // Begin the calculation of B_{i+1}sj
                    while(si_iter!=sj_iter){
                        // Add some additional references to the iterators 
                        const Vector& sj=*sj_iter;
                        Vector& Bisj=*Bisj_iter;

                        // Determine <si,Bisj>
                        Real inner_si_Bisj=VS::innr(si,Bisj);

                        // Determine <yi,sj>
                        Real inner_yi_sj=VS::innr(yi,sj);

                        // Determine -<si,Bisj>/<Bisi,si> Bisi + Bisj
                        // Store in Bisj.  This will become B_{i+1}sj.
                        VS::axpy(-inner_si_Bisj/inner_Bisi_si,Bisi,Bisj);

                        // Determine <yi,sj>/<yi,si> yi + w where we calculated 
                        // w in the line above.  This completes the
                        // computation of B_{i+1}sj.
                        VS::axpy(inner_yi_sj/inner_yi_si,yi,Bisj);

                        // Change j to be j-1 and adjust Bisj and sj accordingly
                        sj_iter++;
                        Bisj_iter++;
                    }

                    // At this point, we've computed all Bisj entries on the
                    // current row.  As a result, we increment i and set j to
                    // be k.  This requires us to modify si, yi, sj, Bisj, and
                    // Bisi accordingly.
                    
                    // Increment i and adjust si
                    si_iter--;

                    // Increment i and adjust yi
                    yi_iter--;

                    // Set j=k and adjust sj
                    sj_iter=oldS.begin();

                    // Set j=k, increment i, and adjust Bisj
                    Bisj_iter=work.begin();

                    // Increment i and adjust Bisi
                    Bisi_iter--;
                }
            }
        };

        // The SR1 Hessian approximation.  
        /* The oldY and oldS lists have the same structure as the BFGS
            preconditioner. */
        class SR1 : public Operator <VS,VS> {
        private:
            const List& oldY;
            const List& oldS;
        public:
            explicit SR1 (const State& state)
                : oldY(state.oldY), oldS(state.oldS) {};
            SR1 (const typename DataStructures <VS>::List& oldY_,
                const typename DataStructures <VS>::List& oldS_)
                : oldY(oldY_), oldS(oldS_) {};
            
            // Operator interface
            void operator () (const Vector& p,Vector& result) const{

                // Check that the number of stored gradient and trial step
                // differences is the same.
                if(oldY.size() != oldS.size())
                    VS::error("In the SR1 Hessian approximation, the number "
                        "of stored gradient differences must equal the number "
                        "of stored trial step differences.");

                // Allocate memory for work
                List work(oldY.size(),p);

                // If we have no vectors in our history, we return the direction
                VS::copy(p,result);
                if(oldY.size() == 0) return;

                // Othwerwise, we copy all of the trial step differences into
                // the work space
                ListIterator Bisj_iter=work.begin();
                ListConstIterator sk_iter=oldS.begin();
                while(Bisj_iter!=work.end())
                    VS::copy((*sk_iter++),(*Bisj_iter++));

                // Keep track of the element Bisi
                ListConstIterator Bisi_iter=work.end(); Bisi_iter--;

                // Keep iterating until Bisi equals the first element in the
                // work list. This means we have computed B1s1, B2s2, ..., Bksk.
                Bisj_iter=work.begin();
                ListConstIterator si_iter=oldS.end(); si_iter--;
                ListConstIterator yi_iter=oldY.end(); yi_iter--;
                ListConstIterator sj_iter=oldS.begin();
                while(1){

                    // Create some reference to our iterators that are easier to
                    // work with
                    const Vector& si=*si_iter;
                    const Vector& yi=*yi_iter;
                    const Vector& Bisi=*Bisi_iter;

                    // Determine <yi,p>
                    Real inner_yi_p=VS::innr(yi,p);

                    // Determine <Bisi,p>
                    Real inner_Bisi_p=VS::innr(Bisi,p);

                    // Determine <yi,si>
                    Real inner_yi_si=VS::innr(yi,si);

                    // Determine <Bisi,si>
                    Real inner_Bisi_si=VS::innr(Bisi,si);

                    // Determine (<yi,p>-<Bisi,p>) / (<y_i,s_i>-<Bisi,si>).
                    // Store in alpha
                    Real alpha=(inner_yi_p-inner_Bisi_p)/(inner_yi_si-inner_Bisi_si);

                    // Determine alpha y_i + Bip.  Store in result (which
                    // accumulate Bip).
                    VS::axpy(alpha,yi,result);

                    // Then, add -alpha*Bisi to this result
                    VS::axpy(-alpha,Bisi,result);

                    // Check whether or not we've calculated B_{i+1}p for the
                    // last time
                    if(Bisi_iter==work.begin()) break;

                    // Begin the calculation of B_{i+1}sj
                    while(si_iter!=sj_iter){
                // Add some additional references to the iterators 
                const Vector& sj=*sj_iter;
                Vector& Bisj=*Bisj_iter;

                // Determine <yi,sj>
                Real inner_yi_sj=VS::innr(yi,sj);

                // Determine <Bisi,sj>
                Real inner_Bisi_sj=VS::innr(Bisi,sj);

                // Determine (<yi,p>-<Bisi,p>) / (<y_i,s_i>-<Bisi,si>).
                // Store in beta 
                Real beta= (inner_yi_sj-inner_Bisi_sj) /
                    (inner_yi_si-inner_Bisi_si);
                
                // Determine beta y_i + Bisj.  Store in Bisj. 
                VS::axpy(beta,yi,Bisj);

                // Add -beta*Bisi to this result
                VS::axpy(-beta,Bisi,Bisj);

                // Change j to be j-1 and adjust Bisj and sj accordingly
                sj_iter++;
                Bisj_iter++;
                    }

                    // At this point, we've computed all Bisj entries on the
                    // current row.  As a result, we increment i and set j to
                    // be k.  This requires us to modify si, yi, sj, Bisj, and
                    // Bisi accordingly.
                    
                    // Increment i and adjust si
                    si_iter--;

                    // Increment i and adjust yi
                    yi_iter--;

                    // Set j=k and adjust sj
                    sj_iter=oldS.begin();

                    // Set j=k, increment i, and adjust Bisj
                    Bisj_iter=work.begin();

                    // Increment i and adjust Bisi
                    Bisi_iter--;
                }
            }
        };

        // The inverse BFGS operator 
        /* The oldY list has the following structure
            oldY[0] = y_k = grad J(u_k) - grad J(u_{k-1})
            oldY[1] = y_{k-1} = grad J(u_{k-1}) - grad J(u_{k-2})
            The oldS list has the following structure
            oldS[0] = s_k = u_k - u_k{-1}
            oldS[1] = s_{k-1} = u_{k-1} - u_k{k-2} */
        class InvBFGS : public Operator <VS,VS> {
        private:
            const List& oldY;
            const List& oldS;
        public:
            explicit InvBFGS(const State& state)
                : oldY(state.oldY), oldS(state.oldS) {};
            
            // Operator interface
            void operator () (const Vector& p,Vector& result) const{

                // Check that the number of stored gradient and trial step
                // differences is the same.
                if(oldY.size() != oldS.size())
                    VS::error("In the inverse BFGS operator, the number "
                        "of stored gradient differences must equal the number "
                        "of stored trial step differences.");
                
                // As a safety check, insure that the inner product between all
                // the (s,y) pairs is positive
                ListConstIterator y0=oldY.begin();
                ListConstIterator s0=oldS.begin();
                while(y0!=oldY.end()){
                    Real inner_y_s=VS::innr(*y0++,*s0++);
                    if(inner_y_s<0)
                        VS::error("Detected a (s,y) pair in the inverse BFGS"
                        "operator that possesed a nonpositive inner product");
                }

                // Create two vectors to hold some intermediate calculations
                std::vector <Real> alpha(oldY.size());
                std::vector <Real> rho(oldY.size());

                // Before we begin computing, copy p to our result 
                VS::copy(p,result);

                // In order to compute, we first iterate over all the stored
                // element in the forward direction.  Then, we iterate over them
                // backward.
                ListConstIterator y_iter=oldY.begin();
                ListConstIterator s_iter=oldS.begin();
                int i=0;
                while(y_iter != oldY.end()){
                    // Find y_k, s_k, and their inner product
                    const Vector& y_k=*(y_iter++);
                    const Vector& s_k=*(s_iter++);
                    rho[i]=Real(1.)/VS::innr(y_k,s_k);

                    // Find rho_i <s_i,result>.  Store in alpha_i
                    alpha[i]=rho[i]*VS::innr(s_k,result);

                    // result = - alpha_i y_i + result 
                    VS::axpy(-alpha[i],y_k,result);

                    // Make sure we don't overwrite alpha and rho
                    i++;
                }

                // Assume that H_0 is the identity operator (which may or may
                // not work in Hilbert space)

                // Now, let us iterate backward over our elements to complete
                // the computation
                while(y_iter != oldY.begin()){
                    // Find y_k and s_k
                    const Vector& s_k=*(--s_iter);
                    const Vector& y_k=*(--y_iter);

                    // beta=rho_i <y_i,result>
                    Real beta= rho[--i] * VS::innr(y_k,result);

                    // result=  (alpha_i-beta) s_i + result
                    VS::axpy(alpha[i]-beta,s_k,result);
                }
            }
        };
        
        // The inverse SR1 operator.  
        /* In this definition, we take a shortcut and simply use the SR1
            Hessian approximation where we swap Y and S.  The oldY and oldS
            lists have the same structure as the BFGS operator. */
        class InvSR1 : public Operator <VS,VS> {
        private:
            SR1 sr1;
        public:
            explicit InvSR1(State& state) : sr1(state.oldS,state.oldY) {};
            void operator () (const Vector& p,Vector& result) const{
                sr1(p,result);
            }
        };

        // Checks a set of stopping conditions
        static StoppingCondition checkStop(const State& state){
            // Create some shortcuts
            const Real& norm_g=state.norm_g;
            const Real& norm_gtyp=state.norm_gtyp;
            const Real& norm_s=state.norm_s;
            const Real& norm_styp=state.norm_styp;
            const int& iter=state.iter;
            const int& iter_max=state.iter_max;
            const Real& eps_g=state.eps_g;
            const Real& eps_s=state.eps_s;

            // Check whether the norm is small relative to some typical gradient
            if(norm_g < eps_g*norm_gtyp)
                return RelativeGradientSmall;

            // Check whether the change in the step length has become too small
            // relative to some typical step
            if(norm_s < eps_s*norm_styp)
                return RelativeStepSmall;

            // Check if we've exceeded the number of iterations
            if(iter>=iter_max)
                return MaxItersExceeded;

            // Otherwise, return that we're not converged 
            return NotConverged;
        }
        
        // Computes the truncated-CG (Steihaug-Toint) trial step for
        // trust-region algorithms
        static void truncatedCG(
            State& state,
            const Operator<VS,VS>& Minv,
            const Operator<VS,VS>& H
        ){

            // Create shortcuts to some elements in the state
            const Vector& u=*(state.u.begin());
            const Vector& g=*(state.g.begin());
            const Real& delta=state.delta;
            const Real& eps_cg=state.eps_krylov;
            const int& iter_max=state.krylov_iter_max;
            Vector& s_k=*(state.s.begin());
            int& iter=state.krylov_iter;
            int& iter_total=state.krylov_iter_total;
            KrylovStop& krylov_stop=state.krylov_stop;
            Real& rel_err=state.krylov_rel_err;

            // Allocate memory for temporaries that we need
            Vector g_k; VS::init(u,g_k);
            Vector v_k; VS::init(u,v_k);
            Vector p_k; VS::init(u,p_k);
            Vector H_pk; VS::init(u,H_pk);

            // Allocate memory for a few constants that we need to track 
            Real kappa;
            Real sigma;
            Real alpha(0.);
            Real beta;
            Real norm_sk_M2,norm_skp1_M2(0.),norm_pk_M2,norm_g;
            Real inner_sk_M_pk,inner_gk_vk,inner_gkp1_vkp1;

            // Initialize our variables
            VS::scal(Real(0.),s_k);                // s_0=0
            VS::copy(g,g_k);                        // g_0=g
            Minv(g_k,v_k);                        // v_0=inv(M)*g_0
            VS::copy(v_k,p_k);                        // p_0=-v_0
            VS::scal(Real(-1.),p_k);
            norm_sk_M2=Real(0.);                // || s_0 ||_M^2 = 0
            norm_pk_M2=VS::innr(g_k,v_k);        // || p_0 ||_M^2 = <g_0,v_0>        
            inner_sk_M_pk=Real(0.);                // <s_0,M p_0>=0
            inner_gk_vk=norm_pk_M2;                // <g_0,v_0> = || p_0 ||_M^2
            norm_g=VS::innr(g,g);                // || g ||

            // Run truncated CG until we hit our max iteration or we converge
            iter_total++;
            for(iter=1;iter<=iter_max;iter++,iter_total++){
                // H_pk=H p_k
                H(p_k,H_pk);

                // Compute the curvature for this direction.  kappa=<p_k,H p_k>
                kappa=VS::innr(p_k,H_pk);

                // If we have negative curvature, don't bother with the next two
                // steps since we're going to exit and we won't need them.  
                if(kappa > 0){
                    // Determine a trial point
                    alpha = VS::innr(g_k,v_k)/kappa;

                    // || s_k+alpha_k p_k ||
                    norm_skp1_M2=norm_sk_M2+Real(2.)*alpha*inner_sk_M_pk
                        +alpha*alpha*norm_pk_M2;
                }

                // If we have negative curvature or our trial point is outside
                // the trust region radius, terminate truncated-CG and find our
                // final step.  We have the kappa!=kappa check in order to trap
                // NaNs.
                if(kappa <= 0 || norm_skp1_M2 >= delta*delta || kappa!=kappa){
                    // sigma = positive root of || s_k + sigma p_k ||_M = delta
                    sigma= (-inner_sk_M_pk + sqrt(inner_sk_M_pk*inner_sk_M_pk
                        + norm_pk_M2*(delta*delta-norm_sk_M2)))/norm_pk_M2;

                    // s_kp1=s_k+sigma p_k
                    VS::axpy(sigma,p_k,s_k);

                    // Return a message as to why we exited
                    if(kappa<=0 || kappa!=kappa)
                        krylov_stop = NegativeCurvature;
                    else
                        krylov_stop = TrustRegionViolated;

                    // Update the residual error for out output,
                    // g_k=g_k+sigma Hp_k
                    VS::axpy(sigma,H_pk,g_k);

                    // Exit the loop
                    break;
                }

                // Take a step in the computed direction. s_k=s_k+alpha p_k
                VS::axpy(alpha,p_k,s_k);

                // Update the norm of sk
                norm_sk_M2=norm_skp1_M2;
                
                // g_k=g_k+alpha H p_k
                VS::axpy(alpha,H_pk,g_k);

                // Test whether we've converged CG
                rel_err=sqrt(VS::innr(g_k,g_k)) / (Real(1e-16)+norm_g);
                if(rel_err <= eps_cg){
                    krylov_stop = RelativeErrorSmall;
                    break;
                }

                // v_k = Minv g_k
                Minv(g_k,v_k);

                // Compute the new <g_kp1,v_kp1>
                inner_gkp1_vkp1=VS::innr(g_k,v_k);

                // beta = <g_kp1,v_kp1> / <g_k,v_k>
                beta= inner_gkp1_vkp1 / inner_gk_vk;

                // Store the new inner product between g_k and p_k
                inner_gk_vk=inner_gkp1_vkp1;
                
                // Find the new search direction.  p_k=-v_k + beta p_k
                VS::scal(beta,p_k);
                VS::axpy(Real(-1.),v_k,p_k);

                // Update the inner product between s_k and M p_k
                inner_sk_M_pk=beta*(inner_sk_M_pk+alpha*norm_pk_M2);

                // Update the norm of p_k
                norm_pk_M2=inner_gk_vk+beta*beta*norm_pk_M2; 

                // Print out diagnostics
                printKrylov(state);
            }

            // Check if we've exceeded the maximum iteration
            if(iter>iter_max){
                krylov_stop=MaxKrylovItersExceeded;
                iter--; iter_total--;
            }
           
            // Grab the relative error in the CG solution
            rel_err=sqrt(VS::innr(g_k,g_k)) / (Real(1e-16)+norm_g);
                
            // Print out diagnostics
            if(iter!=iter_max) printKrylov(state);
        }

        // Checks whether we accept or reject a step
        static bool checkStep(
            State& state,
            const Operator<VS,VS>& H,
            const Functional<VS>& obj_fn
        ){
            // Create shortcuts to some elements in the state
            const Vector& s=*(state.s.begin());
            const Vector& u=*(state.u.begin());
            const Vector& g=*(state.g.begin());
            const Real& eta1=state.eta1;
            const Real& eta2=state.eta2;
            const Real& delta_max=state.delta_max;
            const Real& obj_u=state.obj_u;
            const Real& norm_s=state.norm_s;
            Real& delta=state.delta;
            Real& rho=state.rho;
            Real& obj_ups=state.obj_ups;

            // Allocate memory for temporaries that we need
            Vector ups; VS::init(u,ups);
            Vector Hu_s; VS::init(u,Hu_s);

            // Determine u+s 
            VS::copy(s,ups);
            VS::axpy(Real(1.),u,ups);

            // Determine the objective function evaluated at u+s
            obj_ups=obj_fn(ups);
            
            // Determine H(u)s
            H(s,Hu_s);

            // Determine alpha+<g,s>+.5*<H(u)s,s>
            Real model_s=obj_u+VS::innr(g,s)+Real(.5)*VS::innr(Hu_s,s);

            // Add a safety check in case we don't actually minimize the TR
            // subproblem correctly. This could happen for a variety of reasons.
            // Most notably, if we do not correctly calculate the Hessian
            // approximation, we could have a nonsymmetric approximation.
            // In that case, truncated-CG will exit, but has an undefined
            // result.  In the case that the actual reduction also increases,
            // rho could have an extraneous positive value.  Hence, we require
            // an extra check.
            if(model_s > obj_u){
                delta = norm_s/Real(2.);
                rho = std::numeric_limits<Real>::quiet_NaN(); 
                return false;
            }

            // Determine the ratio of reductions
            rho = (obj_u - obj_ups) / (obj_u - model_s);

            // Update the trust region radius and return whether or not we
            // accept the step
            if(rho >= eta2){
                // Only increase the size of the trust region if we were close
                // to the boundary
                if(fabs(norm_s-delta) < Real(1e-4)*delta)
                    delta = std::min(delta*Real(2.),delta_max);
                return true;
            } else if(rho >= eta1 && rho < eta2)
                return true;
            else {
                delta = norm_s/Real(2.);
                return false;
            }
        }

        // Finds the trust-region step
        static void getStepTR(
            State& state,
            const Operator<VS,VS>& Minv,
            const Operator<VS,VS>& H,
            const Functional<VS>& obj_fn
        ){
            // Create some shortcuts
            const Real& eps_s=state.eps_s;
            int& rejected_trustregion=state.rejected_trustregion;
            Vector& s=*(state.s.begin());
            Real& norm_s=state.norm_s;
            List& oldY=state.oldY;
            List& oldS=state.oldS;
            int& history_reset=state.history_reset;

            // Initialize the counter for the number of rejected steps
            for(rejected_trustregion=0;
                true;
                rejected_trustregion++
            ) {
                // If the number of rejected steps is above the history_reset
                // threshold, destroy the quasi-Newton information
                if(rejected_trustregion > history_reset){
                    oldY.empty();
                    oldS.empty();
                }

                // Print out diagnostic information if we reject a step
                if(rejected_trustregion>0) printState(state,true);

                // Use truncated-CG to find a new trial step
                truncatedCG(state,Minv,H);

                // Save the length of the trial step
                norm_s=sqrt(VS::innr(s,s));
                
                // Check whether the step is good
                if(checkStep(state,H,obj_fn)) break;
                
                // Alternatively, check if the step becomes so small
                // that we're not making progress.  In this case, take
                // a zero step and allow the stopping conditions to exit
                if(norm_s < eps_s) {
                    VS::scal(Real(0.),s);
                    norm_s=Real(0.);
                    break;
                }
            } 
        }

        // Steepest descent search direction
        static void SteepestDescent(State& state){
            // Create some shortcuts 
            const Vector& g=*(state.g.begin());
            Vector& s=*(state.s.begin());

            // We take the steepest descent direction
            VS::copy(g,s);
            VS::scal(Real(-1.),s);
        }

        // Fletcher-Reeves CG search direction
        static void FletcherReeves(State& state){

            // Create some shortcuts 
            const Vector& g=*(state.g.begin());
            const Vector& g_old=*(state.g_old.begin());
            const Vector& s_old=*(state.s_old.begin());
            const int& iter=state.iter;
            Vector& s=*(state.s.begin());

            // If we're on the first iterations, we take the steepest descent
            // direction
            if(iter==1) SteepestDescent(state);

            // On subsequent iterations, we take the FR direction
            else {
                // Find the momentum parameter
                double beta=VS::innr(g,g)/VS::innr(g_old,g_old);

                // Find -g+beta*s_old
                VS::copy(g,s);
                VS::scal(Real(-1.),s);
                VS::axpy(beta,s_old,s);
            }
        }
        
        // Polak-Ribiere CG search direction
        static void PolakRibiere(State& state){

            // Create some shortcuts 
            const Vector& g=*(state.g.begin());
            const Vector& g_old=*(state.g_old.begin());
            const Vector& s_old=*(state.s_old.begin());
            const int& iter=state.iter;
            Vector& s=*(state.s.begin());

            // If we're on the first iterations, we take the steepest descent
            // direction
            if(iter==1) SteepestDescent(state);

            // On subsequent iterations, we take the FR direction
            else {
                // Find the momentum parameter
                double beta=(VS::innr(g,g)-VS::innr(g,g_old))
                    /VS::innr(g_old,g_old);

                // Find -g+beta*s_old
                VS::copy(g,s);
                VS::scal(Real(-1.),s);
                VS::axpy(beta,s_old,s);
            }
        }
        
        // Hestenes-Stiefel search direction
        static void HestenesStiefel(State& state){

            // Create some shortcuts 
            const Vector& g=*(state.g.begin());
            const Vector& g_old=*(state.g_old.begin());
            const Vector& s_old=*(state.s_old.begin());
            const int& iter=state.iter;
            Vector& s=*(state.s.begin());

            // If we're on the first iterations, we take the steepest descent
            // direction
            if(iter==1) SteepestDescent(state);

            // On subsequent iterations, we take the FR direction
            else {
                // Find the momentum parameter
                double beta=(VS::innr(g,g)-VS::innr(g,g_old))
                    /(VS::innr(g,s_old)-VS::innr(g_old,s_old));

                // Find -g+beta*s_old
                VS::copy(g,s);
                VS::scal(Real(-1.),s);
                VS::axpy(beta,s_old,s);
            }
        }

        // BFGS search direction
        static void BFGS_SD(State& state){
            
            // Create some shortcuts 
            const Vector& g=*(state.g.begin());
            Vector& s=*(state.s.begin());

            // Create the inverse BFGS operator
            InvBFGS Hinv(state); 

            // Apply the inverse BFGS operator to the gradient
            Hinv(g,s);

            // Negate the result
            VS::scal(Real(-1.),s);
        }
        
        /* Computes the Newton-CG (truncated-CG) trial step.  Essentially, this
        is the same as trust-region except that we do not have a restriction 
        on the size of the step (no trust-reigon radius).  In the case that we
        encounter negative curvature, we use the last good step.  */ 
        static void NewtonCG(
            State& state,
            const Operator<VS,VS>& Minv,
            const Operator<VS,VS>& H
        ){
            
            // Create shortcuts to some elements in the state
            const Vector& u=*(state.u.begin());
            const Vector& g=*(state.g.begin());
            const Real& eps_cg=state.eps_krylov;
            const int& iter_max=state.krylov_iter_max;
            Vector& s_k=*(state.s.begin());
            int& iter=state.krylov_iter;
            int& iter_total=state.krylov_iter_total;
            KrylovStop& krylov_stop=state.krylov_stop;
            Real& rel_err=state.krylov_rel_err;

            // Allocate memory for temporaries that we need
            Vector g_k; VS::init(u,g_k);
            Vector v_k; VS::init(u,v_k);
            Vector p_k; VS::init(u,p_k);
            Vector H_pk; VS::init(u,H_pk);

            // Allocate memory for a few constants that we need to track 
            Real kappa;
            Real alpha(0.);
            Real beta;
            Real norm_g;
            Real inner_gk_vk,inner_gkp1_vkp1;

            // Initialize our variables
            VS::scal(Real(0.),s_k);                // s_0=0
            VS::copy(g,g_k);                        // g_0=g
            Minv(g_k,v_k);                        // v_0=inv(M)*g_0
            VS::copy(v_k,p_k);                        // p_0=-v_0
            VS::scal(Real(-1.),p_k);
            inner_gk_vk=VS::innr(g_k,v_k);        // <g_0,v_0>        
            norm_g=VS::innr(g,g);                // || g ||

            // Run truncated CG until we hit our max iteration or we converge
            iter_total++;
            for(iter=1;iter<=iter_max;iter++,iter_total++){
                // H_pk=H p_k
                H(p_k,H_pk);

                // Compute the curvature for this direction.  kappa=<p_k,H p_k>
                kappa=VS::innr(p_k,H_pk);

                // If we have negative curvature, don't bother with the next 
                // step since we're going to exit and we don't need it. 
                if(kappa > 0){
                    // Determine a trial point
                    alpha = VS::innr(g_k,v_k)/kappa;
                }

                // If we have negative curvature terminate truncated-CG and find
                // our final step.  We have the kappa!=kappa check in order to
                // trap NaNs.
                if(kappa <= 0 || kappa!=kappa){

                    // If we're on the first iteration and we already have
                    // negative curvature, use the steepest-descent direction.
                    if(iter==1){
                        VS::copy(g_k,s_k);
                        VS::scal(Real(-1.),s_k);
                    }

                    // Return a message as to why we exited
                    krylov_stop = NegativeCurvature;

                    // Exit the loop
                    break;
                }

                // Take a step in the computed direction. s_k=s_k+alpha p_k
                VS::axpy(alpha,p_k,s_k);

                // g_k=g_k+alpha H p_k
                VS::axpy(alpha,H_pk,g_k);

                // Test whether we've converged CG
                if(sqrt(VS::innr(g_k,g_k)) <= eps_cg*norm_g){
                    krylov_stop = RelativeErrorSmall;
                    break;
                }

                // v_k = Minv g_k
                Minv(g_k,v_k);

                // Compute the new <g_kp1,v_kp1>
                inner_gkp1_vkp1=VS::innr(g_k,v_k);

                // beta = <g_kp1,v_kp1> / <g_k,v_k>
                beta= inner_gkp1_vkp1 / inner_gk_vk;

                // Store the new inner product between g_k and p_k
                inner_gk_vk=inner_gkp1_vkp1;
                
                // Find the new search direction.  p_k=-v_k + beta p_k
                VS::scal(beta,p_k);
                VS::axpy(Real(-1.),v_k,p_k);
            }

            // Check if we've exceeded the maximum iteration
            if(iter>iter_max){
              krylov_stop=MaxKrylovItersExceeded;
              iter--; iter_total--;
            }
           
            // Grab the relative error in the CG solution
            rel_err=sqrt(VS::innr(g_k,g_k)) / (Real(1e-16)+norm_g);
        }

        // Compute a Golden-Section search between eps and 2*alpha where
        // alpha is the last line search parameter.
        static void goldenSection(
            State& state,
            const Functional<VS>& F
        ) {
            // Create some shortcuts
            const Vector& u=*(state.u.begin());
            const int& iter_max=state.linesearch_iter_max;
            Real& alpha=state.alpha;
            Vector& s=*(state.s.begin());
            int& iter_total=state.linesearch_iter_total;
            int& iter=state.linesearch_iter;
            Real& obj_ups=state.obj_ups;

            // Create one work element
            Vector work; VS::init(u,work);

            // Find 1 over the golden ratio
            Real beta=Real(2./(1.+sqrt(5.)));

            // Find a bracket for the linesearch such that a < b
            Real a=Real(1e-16);
            Real b=Real(2.)*alpha;

            // Find two new points between a and b, mu and lambda,
            // such that lambda < mu
            double lambda=a+(1.-beta)*(b-a);
            double mu=a+beta*(b-a);

            // Find the objective value at mu and labmda 

            // mu 
            VS::copy(u,work);
            VS::axpy(mu,s,work);
            Real obj_mu=F(work);

            // lambda
            VS::copy(u,work);
            VS::axpy(lambda,s,work);
            Real obj_lambda=F(work);

            // Search for a fixed number of iterations 
            for(iter=1;iter<=iter_max;iter++,iter_total++){

                // If the objective is greater on the left, bracket on the
                // right.  Alternatively, it's possible that we're going to
                // generate a NaN on the right.  This means that obj_mu=NaN.
                // In this case we want to bracket on the left.  Since
                // obj_lambda > obj_mu will return false when obj_mu is a NaN,
                // we should be safe.
                if(obj_lambda > obj_mu){
                    a=lambda;
                    lambda=mu;
                    obj_lambda=obj_mu;
                    mu=a+beta*(b-a);

                    VS::copy(u,work);
                    VS::axpy(mu,s,work);
                    obj_mu=F(work);

                // Otherwise, the objective is greater on the right, so bracket
                // on the left
                } else {
                    b=mu;
                    mu=lambda;
                    obj_mu=obj_lambda;
                    lambda=a+(1-beta)*(b-a);
            
                    VS::copy(u,work);
                    VS::axpy(lambda,s,work);
                    obj_lambda=F(work);
                }
            }

            // The iteration count is technically one larger than it should be
            iter--; iter_total--;

            // Once we're finished narrowing in on a solution, take our best
            // guess for the line search parameter
            alpha=obj_lambda < obj_mu ? lambda : mu;

            // Save the objective value at this step
            obj_ups=obj_lambda < obj_mu ? obj_lambda : obj_mu;
        }

        // Find the line search parameter based on the 2-point approximation
        // from Barzilai and Borwein
        static void twoPoint(
            State& state,
            const Functional<VS>& F
        ) {
            // Create some shortcuts
            const Vector& u=*(state.u.begin());
            const Vector& g=*(state.g.begin());
            const Vector& u_old=*(state.u_old.begin());
            const Vector& g_old=*(state.g_old.begin());
            const LineSearchKind& kind=state.kind;
            Real& alpha=state.alpha;
            Vector& s=*(state.s.begin());
            int& iter_total=state.linesearch_iter_total;
            int& iter=state.linesearch_iter;
            Real& obj_ups=state.obj_ups;

            // Create elements for delta_u and delta_g as well as one work
            // element
            Vector delta_u; VS::init(u,delta_u);
            Vector delta_g; VS::init(u,delta_g);
            Vector work; VS::init(u,work);

            // Find delta_u
            VS::copy(u,delta_u);
            VS::axpy(Real(-1.),u_old,delta_u);

            // Find delta_g
            VS::copy(g,delta_g);
            VS::axpy(Real(-1.),g_old,delta_g);

            // Find alpha
            if(kind==TwoPointA_t)
                alpha=VS::innr(delta_u,delta_g)/VS::innr(delta_g,delta_g);
            else if(kind==TwoPointB_t)
                alpha=VS::innr(delta_u,delta_u)/VS::innr(delta_u,delta_g);

            // Save the objective value at this step
            VS::copy(u,work);
            VS::axpy(alpha,s,work);
            obj_ups=F(work);

            // Since we do one function evaluation, increase the linesearch
            // iteration by one
            iter=1; iter_total++;
        }
        
        // Compute a backtracking line-search. 
        static void backTracking(
            State& state,
            const Functional<VS>& F
        ) {
            // Create some shortcuts
            const Vector& u=*(state.u.begin());
            const int& iter_max=state.linesearch_iter_max;
            Real& alpha=state.alpha;
            Vector& s=*(state.s.begin());
            int& iter_total=state.linesearch_iter_total;
            int& iter=state.linesearch_iter;
            Real& obj_ups=state.obj_ups;

            // Create one work element
            Vector work; VS::init(u,work);

            // Store the best objective value and alpha that we used to find it.
            // Our initial guess will be at alpha*2.
            Real alpha_best=Real(2.)*alpha;
            VS::copy(u,work);
            VS::axpy(alpha_best,s,work);
            Real obj_best=F(work);

            // Evaluate the objective iter_max times at a distance of
            // 2*alpha, alpha, alpha/2, ....  Then, pick the best one.
            Real alpha0=alpha;
            for(int i=0;i<iter_max-1;i++){
                    // Evaluate F(u+alpha*s)
                VS::copy(u,work);
                VS::axpy(alpha0,s,work);
                Real obj=F(work);

                // If this is better than our best guess so far, save it
                if(obj<obj_best){
                    obj_best=obj;
                    alpha_best=alpha0;
                }

                // Reduce the size of alpha
                alpha0 /= Real(2.);
            }

            // Save the best objective and alpha found
            alpha=alpha_best;
            obj_ups=obj_best;

            // Indicate how many iterations we used to find this value
            iter_total+=iter_max;
            iter=iter_max;
        }
       
        // Finds a trial step using a line-search for globalization
        static void getStepLS(
            State& state,
            const Operator<VS,VS>& Minv,
            const Operator<VS,VS>& H,
            const Functional<VS>& obj_fn
        ){
            // Create some shortcuts
            const LineSearchDirection& dir=state.dir;
            const LineSearchKind& kind=state.kind;
            const int& iter=state.iter;
            const int& linesearch_iter_max=state.linesearch_iter_max;
            const Real& obj_u=state.obj_u;
            const Real& obj_ups=state.obj_ups;
            const Real& eps_s=state.eps_s;
            Vector& s=*(state.s.begin());
            Real& norm_s=state.norm_s;
            Real& alpha=state.alpha;

            // Find the line-search direction
            switch(dir){
            case SteepestDescent_t:
                SteepestDescent(state);
                break;
            case FletcherReeves_t:
                FletcherReeves(state);
                break;
            case PolakRibiere_t:
                PolakRibiere(state);
                break;
            case HestenesStiefel_t:
                HestenesStiefel(state);
                break;
            case LimitedMemoryBFGS_t:
                BFGS_SD(state);
                break;
            case NewtonCG_t:
                NewtonCG(state,Minv,H);
                break;
            }

            // Do a line-search in the specified direction
            switch(kind){
            case GoldenSection_t:
                do{
                    // Conduct the golden section search
                    goldenSection(state,obj_fn);

                    // If we don't decrease, print out some diagnostic
                    // information and reduce the size of alpha
                    if(obj_ups > obj_u){
                        norm_s=alpha*sqrt(VS::innr(s,s));
                        printState(state,true);
                        alpha /= Real(4.);

                        // Check if the step becomes so small that we're not
                        // making progress.  In this case, take a zero step 
                        // and allow the stopping conditions to exit
                        if(norm_s < eps_s) {
                            alpha=0.;
                            break;
                        }
                    }

                // If we don't decrease the objective, try another linesearch
                } while(obj_u < obj_ups || obj_ups!=obj_ups);
                break;
            case BackTracking_t:
                do{
                    // Conduct a backtracking search
                    backTracking(state,obj_fn);

                    // If we don't decrease, print out some diagnostic
                    // information and restart the search at the smallest
                    // alpha we previously searched.
                    if(obj_ups > obj_u){
                        norm_s=alpha*sqrt(VS::innr(s,s));
                        printState(state,true);
                        alpha = alpha/pow(Real(2.),linesearch_iter_max+1);

                        // Check if the step becomes so small that we're not
                        // making progress.  In this case, take a zero step 
                        // and allow the stopping conditions to exit
                        if(norm_s < eps_s) {
                            alpha=0.;
                            break;
                        }
                    }

                // If we don't decrease the objective, try another linesearch
                } while(obj_u < obj_ups || obj_ups!=obj_ups);
                break;
            case TwoPointA_t:
            case TwoPointB_t:
                    if(iter>1) twoPoint(state,obj_fn);
                else goldenSection(state,obj_fn);
                break;
            case Brents_t:
                VS::error("Brent's linesearch not currently implemented.");
            }
        
            // Scale the line-search direction by the line search parameter 
            VS::scal(alpha,s);

            // Save the length of the trial step
            norm_s=sqrt(VS::innr(s,s));
        }

        // Finds a new trial step
        static void getStep(
            State& state,
            const Operator<VS,VS>& Minv,
            const Operator<VS,VS>& H,
            const Functional<VS>& obj_fn
        ){
            // Create some shortcuts
            const AlgorithmClass& algorithm_class=state.algorithm_class;

            // Choose whether we use a line-search or trust-region method
            switch(algorithm_class){
            case TrustRegion:
                getStepTR(state,Minv,H,obj_fn);
                break;
            case LineSearch:
                getStepLS(state,Minv,H,obj_fn);
                break;
            }
        }

        // Updates the quasi-Newton information
        static void updateQuasi(State& state) {
            // Exit immediately if we're not using a quasi-Newton method
            if(state.stored_history==0) return;

            // Create some shortcuts
            const Vector& u=*(state.u.begin());
            const Vector& g=*(state.g.begin());
            const Vector& u_old=*(state.u_old.begin());
            const Vector& g_old=*(state.g_old.begin());
            const Operators& Minv_type=state.Minv_type;
            const Operators& H_type=state.H_type;
            const LineSearchDirection& dir=state.dir;
            List& oldY=state.oldY;
            List& oldS=state.oldS;
           
            // Allocate some temp storage for y and s
            Vector s; VS::init(u,s);
            Vector y; VS::init(u,y);

            // Find s = u-u_old
            VS::copy(u,s);
            VS::axpy(Real(-1.),u_old,s);

            // Find y = g - g_old
            VS::copy(g,y);
            VS::axpy(Real(-1.),g_old,y);

            // If we're using BFGS, check that <y,x> > 0
            if((Minv_type==InvBFGS_t || H_type==BFGS_t
                || dir==LimitedMemoryBFGS_t)
                && VS::innr(y,s) < 0)
                return;

            // Insert these into the quasi-Newton storage
            oldS.push_front(s);
            oldY.push_front(y);

            // Determine if we need to free some memory
            if(oldS.size()>state.stored_history){
                    oldS.pop_back();
                oldY.pop_back();
            }
        }

        // Solves an optimization problem
        static void getMin(
            State& state,
            StateManipulator& smanip,
            Functional<VS>& F,
            Operator<VS,VS>& G,
            Operator<VS,VS>& H,
            Operator<VS,VS>& Minv
        ) {
            // Create some shortcuts
            Vector& u=*(state.u.begin());
            Vector& g=*(state.g.begin());
            Vector& s=*(state.s.begin());
            Vector& u_old=*(state.u_old.begin());
            Vector& g_old=*(state.g_old.begin());
            Vector& s_old=*(state.s_old.begin());
            Real& obj_u=state.obj_u;
            Real& obj_ups=state.obj_ups;
            Real& norm_s=state.norm_s;
            Real& norm_g=state.norm_g;
            Real& norm_gtyp=state.norm_gtyp;
            Real& norm_styp=state.norm_styp;
            int& iter=state.iter;
            StoppingCondition& opt_stop=state.opt_stop;
            AlgorithmClass& algorithm_class=state.algorithm_class;
            LineSearchDirection& dir=state.dir;

                    
            // Evaluate the objective function and gradient
            obj_u=F(u);
            G(u,g);
            norm_g=sqrt(VS::innr(g,g));
            norm_gtyp=norm_g;

            // Prints out the header for the diagonstic information
            printStateHeader(state);
            if(algorithm_class==TrustRegion || dir==NewtonCG_t)
                printKrylovHeader(state);

            // Primary optimization loop
            do{
                    // Print some diagnostic information
                printState(state);

                // Get a new optimization iterate.  
                getStep(state,Minv,H,F);

                // On the first iteration, save the size of the step
                if(iter==1) norm_styp=norm_s;

                // Save the old variable, gradient, and trial step.  This
                // is useful for both CG and quasi-Newton methods.
                VS::copy(u,u_old);
                VS::copy(g,g_old);
                VS::copy(s,s_old);

                // Move to the new iterate
                VS::axpy(Real(1.),s,u);

                // Manipulate the state if required
                smanip(state);

                // Find the new objective function and gradient
                obj_u=obj_ups;
                G(u,g);
                norm_g=sqrt(VS::innr(g,g));

                // Update the quasi-Newton information
                updateQuasi(state);

                // Increase the iteration
                iter++;
            } while((opt_stop=checkStop(state))==NotConverged);
                    
            // Print a final diagnostic 
            printState(state);
        }
        
        // Solves an optimization problem
        static void getMin(
            State& state,
            StateManipulator& smanip,
            Functional<VS>& F,
            Operator<VS,VS>& G,
            Operator<VS,VS>& H
        ){
            // Create a few shortcuts
            Operators& Minv_type=state.Minv_type;

            // Determine the preconditioner
            std::auto_ptr <Operator<VS,VS> > Minv;
            switch(Minv_type){
                case Identity_t:
                    Minv.reset(new Identity());
                    break;
                case InvBFGS_t:
                    Minv.reset(new InvBFGS(state));
                    break;
                case InvSR1_t:
                    Minv.reset(new InvSR1(state));
                    break;
                case External_t:
                    VS::error("An externally defined preconditioner must be "
                        "provided explicitely.");
                    break;
                default:
                    VS::error("Not a valid Hessian approximation.");
                    break;
            }

            // Solve the minimization problem
            getMin(state,smanip,F,G,H,*Minv);
        }
        
        // Solves an optimization problem
        static void getMin(
            State& state,
            Functional<VS>& F,
            Operator<VS,VS>& G,
            Operator<VS,VS>& H,
            Operator<VS,VS>& Minv
        ) {
            StateManipulator smanip;
            getMin(state,smanip,F,G,H,Minv);
        }
        
        // Solves an optimization problem
        static void getMin(
            State& state,
            Functional<VS>& F,
            Operator<VS,VS>& G,
            Operator<VS,VS>& H
        ){
            StateManipulator smanip;
            getMin(state,smanip,F,G,H);
        }
        
        // Solves an optimization problem
        static void getMin(
            State& state,
            StateManipulator& smanip,
            Functional<VS>& F,
            Operator<VS,VS>& G
        ){
            // Create a few shortcuts
            Operators& H_type=state.H_type;

            // Determine the Hessian approximation
            std::auto_ptr <Operator<VS,VS> > H;
            switch(H_type){
                    case Identity_t:
                    H.reset(new Identity());
                    break;
                    case ScaledIdentity_t:
                    H.reset(new ScaledIdentity(state));
                    break;
                case BFGS_t:
                    H.reset(new BFGS(state));
                    break;
                case SR1_t:
                    H.reset(new SR1(state));
                    break;
                case External_t:
                    VS::error("An externally defined Hessian must be provided "
                        "explicitely.");
                    break;
                default:
                    VS::error("Not a valid Hessian approximation.");
                    break;
            }

            // Solve the minimization problem
            getMin(state,smanip,F,G,*H);
        }
        
        // Solves an optimization problem
        static void getMin(
            State& state,
            Functional<VS>& F,
            Operator<VS,VS>& G
        ){
            StateManipulator smanip;
            getMin(state,smanip,F,G);
        }

        // Prints out useful information regarding the current optimization
        // state
        static void printState(const State& state,const bool noiter=false){
            // Create some shortcuts
            const int& iter=state.iter;
            const Real& obj_u=state.obj_u;
            const Real& norm_g=state.norm_g;
            const Real& norm_s=state.norm_s;
            const Real& krylov_rel_err=state.krylov_rel_err;
            const int& krylov_iter=state.krylov_iter;
            const KrylovStop& krylov_stop=state.krylov_stop; 
            const AlgorithmClass& algorithm_class=state.algorithm_class;
            const LineSearchDirection& dir=state.dir;
            const int& linesearch_iter=state.linesearch_iter;
            const int& verbose=state.verbose;

            // Check if we should print
            if(verbose<1) return;

            // Basic information
            std::stringstream ss;
            if(!noiter) ss << std::setw(4) << iter << ' ';
            else ss << std::setw(4) << '*' << ' ';
            ss << std::scientific << std::setprecision(3)
                    << std::setw(11) << obj_u << ' '
                << std::setw(11) << norm_g << ' ';
            if(iter==0) ss << "            ";
            else ss << std::setw(11) << norm_s << ' ';

            // Information for the Krylov method
            if(algorithm_class==TrustRegion || dir==NewtonCG_t){
                ss << std::setw(11) << krylov_rel_err << ' ' 
                    << std::setw(6) << krylov_iter << ' ';

                switch(krylov_stop){
                case NegativeCurvature:
                    ss << std::setw(10) << "Neg Curv" << ' '; 
                    break;
                case RelativeErrorSmall:
                    ss << std::setw(10) << "Rel Err " << ' '; 
                    break;
                case MaxKrylovItersExceeded:
                    ss << std::setw(10) << "Max Iter" << ' '; 
                    break;
                case TrustRegionViolated:
                    ss << std::setw(10) << "Trst Reg" << ' '; 
                    break;
                }
            }

            // Information for the line search
            if(algorithm_class==LineSearch){
                ss << std::setw(6) << linesearch_iter << ' ';
            }

            // Send the information to the screen
            ss << std::endl;
            VS::print(ss.str());
        }
        
        // Prints out useful information regarding the Krylov method 
        static void printKrylov(const State& state){
            // Create some shortcuts
            const Real& krylov_rel_err=state.krylov_rel_err;
            const int& krylov_iter=state.krylov_iter;
            const int& krylov_iter_total=state.krylov_iter_total;
            const int& verbose=state.verbose;

            // Check if we should print
            if(verbose<2) return;

            // Basic information
            std::stringstream ss;
            ss << "  " << std::setw(4) << krylov_iter << ' '
                << std::setw(6) << krylov_iter_total << ' '
                << std::scientific << std::setprecision(3)
                    << std::setw(11) << krylov_rel_err; 

            // Send the information to the screen
            ss << std::endl;
            VS::print(ss.str());
        }

        // Prints out a description for the state
        static void printStateHeader(const State& state){
            // Create some shortcuts
            const AlgorithmClass& algorithm_class=state.algorithm_class;
            const LineSearchDirection& dir=state.dir;
            const int& verbose=state.verbose;

            // Check if we should print
            if(verbose<1) return;

            // Basic Header
            std::stringstream ss;
            ss << "Iter" << ' '
                    << std::setw(11) << "Obj Value" << ' '
                    << std::setw(11) << "norm(g)  " << ' '
                    << std::setw(11) << "norm(s)  " << ' ';

            // In case we're using a Krylov method
            if(algorithm_class==TrustRegion || dir==NewtonCG_t){
                    ss << std::setw(11) << "Kry Error" << ' '
                    << std::setw(6) << "KryIt" << ' ' 
                    << std::setw(10) << "Kry Why " << ' ';
            }

            // Information for the line search
            if(algorithm_class==LineSearch){
                ss << std::setw(6) << "LS It" << ' ';
            }

            // Send the information to the screen
            ss << std::endl;
            VS::print(ss.str());
        }
        
        // Prints out a description for the Krylov method 
        static void printKrylovHeader(const State& state){
            // Create some shortcuts
            const int& verbose=state.verbose;

            // Check if we should print
            if(verbose<2) return;

            // Basic Header
            std::stringstream ss;
            ss << "  Iter" << ' '
                    << std::setw(6) << "Total" << ' '
                    << std::setw(11) << "Rel Err" << ' ';

            // Send the information to the screen
            ss << std::endl;
            VS::print(ss.str());
        }
    };


    // Utilities for parameter estimation.  
    template <class U,class Y,class Z>
    class parest{
    private:
        // Prevent instantiation of this class
        explicit parest();
        // No assignment operator required
        // No copy constructor required
    public:
        // Setup some types
        typedef typename U::Vector U_Vector; 
        typedef typename U::Real U_Real; 
        typedef typename Y::Vector Y_Vector; 
        typedef typename Y::Real Y_Real; 
        typedef typename Z::Vector Z_Vector; 

        // The Gauss-Newton Hessian approximation
        class GaussNewton : public Operator <U,U> {
        private:
            // Residual function
            const DiffOperator <Y,Z>& r;

            // Solution operator
            const DiffOperator <U,Y>& h;

            // Point around which we find the Hessian approximation
            const U_Vector& u;

            // Points in the Y and Z space that we use to initialize
            // additional memory.
            const Y_Vector& y;
            const Z_Vector& z;
        public:
            GaussNewton(
                const DiffOperator <Y,Z>& r_,
                const DiffOperator <U,Y>& h_,
                const U_Vector& u_,
                const Y_Vector& y_,
                const Z_Vector& z_) :
                r(r_), h(h_), u(u_), y(y_), z(z_)
            {
                    // Insure the max index on r and h are the same
                if(r.max_index()!=h.max_index())
                    U::error("The solution and residual operators used "
                            "to compose the Gauss-Newton operator have "
                        "a differing max index.");
            }
           
            // Hessian vector product
            void operator () (
                const U_Vector& p,
                U_Vector& result
            ) const {

                // Create two work elements in Y, one in Z, and one in U
                Y_Vector y1; Y::init(y,y1);
                Y_Vector y2; Y::init(y,y2);
                Z_Vector z1; Z::init(z,z1);
                U_Vector u1; U::init(u,u1);

                // Zero out the result 
                U::scal(U_Real(0.),result);

                // Accumulate the affect of the GN approximation on p one piece
                // at a time
                for(int i=0;i<r.max_index();i++){

                    // Find the solution, y.  Store in y1. 
                    (*h.f(i))(u,y1);

                    // Find h'(u)p.  Store in y2. 
                    (*h.fp(i,u))(p,y2);

                    // Find r'(y) y2.  Store in z1.
                    (*r.fp(i,y1))(y2,z1);        

                    // Find f'(y)* z1.  Store in y2. 
                    (*r.fps(i,y1))(z1,y2);

                    // Find h'(u)* y2.  Store in u1. 
                    (*h.fps(i,u))(y2,u1);

                    // Accumulate the result
                    U::axpy(1.,u1,result);
                }
            }
        };
        
        // The full-Newton Hessian 
        class Newton : public Operator <U,U> {
        private:
            // Residual function
            const DiffOperator <Y,Z>& r;

            // Solution operator
            const DiffOperator <U,Y>& h;

            // Point around which we find the Hessian approximation
            const U_Vector& u;

            // Points in the Y and Z space that we use to initialize
            // additional memory.
            const Y_Vector& y;
            const Z_Vector& z;

        public:
            Newton(
                const DiffOperator <Y,Z>& r_,
                const DiffOperator <U,Y>& h_,
                const U_Vector& u_,
                const Y_Vector& y_,
                const Z_Vector& z_) :
                r(r_), h(h_), u(u_), y(y_), z(z_)
            {
                    // Insure the max index on r and h are the same
                if(r.max_index()!=h.max_index())
                    U::error("The solution and residual operators used "
                            "to compose the Newton operator have "
                        "a differing max index.");
            }

            // Operator-vector product
            void operator () (
                const U_Vector& p,
                U_Vector& result
            ) const{

                // Create three work elements in Y, one in Z, and one in U
                Y_Vector y1; Y::init(y,y1);
                Y_Vector y2; Y::init(y,y2);
                Y_Vector y3; Y::init(y,y3);
                Z_Vector z1; Z::init(z,z1);
                U_Vector u1; U::init(u,u1);

                // Zero out the result 
                U::scal(U_Real(0.),result);

                // Accumulate the affect of the GN approximation on p one piece
                // at a time
                for(int i=0;i<r.max_index();i++){

                    // Accumulate the Gauss-Newton part of the Hessian

                    // Find the solution, y.  Store in y1.  SAVE THIS. 
                    (*h.f(i))(u,y1);

                    // Find h'(u)p.  Store in y2.  SAVE THIS. 
                    (*h.fp(i,u))(p,y2);

                    // Find r'(y) y2.  Store in z1.
                    (*r.fp(i,y1))(y2,z1);        

                    // Find f'(y)* z1.  Store in y3. 
                    (*r.fps(i,y1))(z1,y3);

                    // Find h'(u)* y2.  Store in u1. 
                    (*h.fps(i,u))(y3,u1);

                    // Accumulate the result
                    U::axpy(U_Real(1.),u1,result);

                    // At this point, y1=h(u) and y2=h'(u)p.
                    
                    // Accumlate the second order terms for the solution
                    // operator

                    // Find r(y).  Store in z1.  SAVE THIS.
                    (*r.f(i))(y1,z1);

                    // Find r'(y)* z1.  Store in y3. 
                    (*r.fps(i,y1))(z1,y3);

                    // Find (h''(u)p)*y3.  Store in u1.
                    (*h.fpps(i,u,y3))(p,u1);

                    // Accumulate the result
                    U::axpy(U_Real(1.),u1,result);

                    // Accumulate the second order terms for the
                    // projection operator

                    // At this point, y1=h(u), y2=h'(u)p, z1=r(y).

                    // Find (r''(y) y2)* z1.  Store in y3.
                    (*r.fpps(i,y1,z1))(y2,y3);

                    // Find h'(u)* y3.  Store in u1. 
                    (*h.fps(i,u))(y3,u1);

                    // Accumulate the result
                    U::axpy(U_Real(1.),u1,result);
                }
            }
        };
        
        // Finds the gradient 
        class getGradient: public Operator <U,U> {
        private:
            // Residual function
            const DiffOperator <Y,Z>& r;

            // Solution operator
            const DiffOperator <U,Y>& h;

            // Points in the Y and Z space that we use to initialize
            // additional memory.
            const Y_Vector& y;
            const Z_Vector& z;

        public:
            getGradient(
                const DiffOperator <Y,Z>& r_,
                const DiffOperator <U,Y>& h_,
                const Y_Vector& y_,
                const Z_Vector& z_) :
                r(r_), h(h_), y(y_), z(z_)
            {
                    // Insure the max index on r and h are the same
                if(r.max_index()!=h.max_index())
                    U::error("The solution and residual operators used "
                            "to compose the gradient operator have "
                        "a differing max index.");
            }
            
            void operator () (const U_Vector& u,U_Vector& g) const {
                // Create two work elements in Y, one in Z, and one in U
                Y_Vector y1; Y::init(y,y1);
                Y_Vector y2; Y::init(y,y2);
                Z_Vector z1; Z::init(z,z1);
                U_Vector u1; U::init(u,u1);

                // Zero out the gradient
                U::scal(U_Real(0.),g);

                // Accumulate the gradient one piece at a time
                for(int i=0;i<r.max_index();i++){

                    // Find the solution y.  Store in y1 
                    (*h.f(i))(u,y1);

                    // Find r(y).  Store in z1.
                    (*r.f(i))(y1,z1);

                    // Find r'(y)* z1.  Store in y2. 
                    (*r.fps(i,y1))(z1,y2);

                    // Find h'(u)* y2.  Store in u1. 
                    (*h.fps(i,u))(y2,u1);

                    // Accumulate the result 
                    U::axpy(U_Real(1.),u1,g);
                }
            }
        };

        // Computes the objective value 
        class getObjective : public Functional <U> {
        private:
            // Residual function
            const DiffOperator <Y,Z>& r;

            // Solution operator
            const DiffOperator <U,Y>& h;

            // Points in the Y and Z space that we use to initialize
            // additional memory.
            const Y_Vector& y;
            const Z_Vector& z;

        public:
            getObjective(
                const DiffOperator <Y,Z>& r_,
                const DiffOperator <U,Y>& h_,
                const Y_Vector& y_,
                const Z_Vector& z_) :
                r(r_), h(h_), y(y_), z(z_)
            {
                    // Insure the max index on r and h are the same
                if(r.max_index()!=h.max_index())
                    U::error("The solution and residual operators used "
                            "to compose the objective function have "
                        "a differing max index.");
            }

            // Finds the objective value
            typename U::Real operator () (const U_Vector& u) const{
                // Create one work element in Y and one in Z
                Y_Vector y1; Y::init(y,y1);
                Z_Vector z1; Z::init(z,z1);

                // Accumulate the objective value a piece at a time
                double obj_u=0.;
                for(int i=0; i<r.max_index(); i++){

                    // Find the solution y.  Store in y1.
                    (*h.f(i))(u,y1);

                    // Find r(y).  Store it in z1. 
                    (*r.f(i))(y1,z1);

                    // Accumulate .5*|| r_i(h_i(u)) ||^2
                    obj_u += .5*Z::innr(z1,z1);
                }

                // Return the accumulated answer
                return obj_u;
            }
        };
    };
}
#endif
