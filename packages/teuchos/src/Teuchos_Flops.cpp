// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Flops.hpp"

namespace Teuchos
{

Flops::Flops(void) : flops_(0.0)
{
}

Flops::Flops(const Flops& flops) : flops_(0.0)
{
}

Flops::~Flops(void)  
{
  flops_ = 0.0;
}

}  // namespace Teuchos
