#include <stk_mms/stk_mms.h>
#include <stk_util/diag/UserPlugin.hpp>

namespace stk
{
namespace mms
{
  // global params
  const double eps = 0.05;

  const double gamma = 1.4;
  const double Rgas = 287.0973757;
  const double Cv = Rgas/(gamma-1.0);
  
  const double M0 = 2.5;
  const double p0 = 35651.28116;
  const double T0 = 236.215;
  
  const double rho0 = p0/(Rgas*T0);
  const double u0 = M0*sqrt(gamma*p0/rho0);

extern "C" FAD2_Type euler_p(const FAD2_Type & x,
			    const FAD2_Type & y,
			    const FAD2_Type & z,
			    const FAD2_Type & t) {
  return p0*(1.0 - eps*sin(M_PI*x)*cos(M_PI*y));
}

extern "C" FAD2_Type euler_T(const FAD2_Type & x,
			    const FAD2_Type & y,
			    const FAD2_Type & z,
			    const FAD2_Type & t) {
  return T0*(1.0 + eps*sin(M_PI*x)*sin(M_PI*y));
}

extern "C" FAD2_Type euler_u(const FAD2_Type & x,
			    const FAD2_Type & y,
			    const FAD2_Type & z,
			    const FAD2_Type & t) {
  return u0*(1.0 - eps*sin(M_PI*x)*cos(M_PI*y));
}

extern "C" FAD2_Type euler_v(const FAD2_Type & x,
			    const FAD2_Type & y,
			    const FAD2_Type & z,
			    const FAD2_Type & t) {
  return u0*eps*sin(M_PI*x)*sin(M_PI*y);
}

extern "C" void register_euler()
{
  sierra::Plugin::UserSubroutine<scalar_FAD2_function>::instance().
    registerFunction(std::string("euler_p"), euler_p);

  sierra::Plugin::UserSubroutine<scalar_FAD2_function>::instance().
    registerFunction(std::string("euler_T"), euler_T);

  sierra::Plugin::UserSubroutine<scalar_FAD2_function>::instance().
    registerFunction(std::string("euler_u"), euler_u);

  sierra::Plugin::UserSubroutine<scalar_FAD2_function>::instance().
    registerFunction(std::string("euler_v"), euler_v);
}

}
}
