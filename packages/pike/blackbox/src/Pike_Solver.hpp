#ifndef PIKE_SOLVER_HPP
#define PIKE_SOLVER_HPP

namespace pike {

  class Solver {

    virtual ~Solver() {}

    virtual void step() = 0;

    virtual void solve() = 0;

  };

}

#endif
