// Copyright (c) 2014-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include "apr_util.h"     // for conv_string
#include "aprepro.h"      // for Aprepro, symrec, etc
#include "init_structs.h" // for svar_init, var_init
#include <cstddef>        // for size_t
#include <cstring>        // for strcmp
#include <iomanip>        // for operator<<, setw, etc
#include <ostream>        // for operator<<, basic_ostream, etc
#include <string>         // for char_traits, operator<<, etc

namespace SEAMS {

  extern Aprepro *aprepro;
  extern int      echo;

  const char *do_Units(char *type);
  void load_conversion(struct var_init *base, struct svar_init *label);

#define DEFINE_VAR(name, val, label)                                                               \
  do {                                                                                             \
    if ((ptr = aprepro->getsym((name))) == nullptr)                                                \
      ptr          = aprepro->putsym((name), SEAMS::Aprepro::VARIABLE, 1);                         \
    ptr->value.var = (val);                                                                        \
    if (echo) {                                                                                    \
      *(aprepro->infoStream) << comment << " 1 " << std::left << std::setw(10) << name             \
                             << "\t= " << std::setw(14) << std::setprecision(7) << val << "  "     \
                             << label << '\n';                                                     \
    }                                                                                              \
  } while (0)

  // clang-format off
namespace {
/*-------------------------------------------------------------------------------------*/
/* SI Units */
struct svar_init si_label[] =
  {
    {"tout", "second"},
    {"lout", "meter"},
    {"aout", "m/sec^2"},
    {"mout", "kilogram"},
    {"fout", "newton"},
    {"vout", "meter/sec"},
    {"Vout", "meter^3"},
    {"dout", "kg/m^3"},
    {"eout", "joule (Nm)"},
    {"Pout", "watt (Nm/s)"},
    {"pout", "Pa"},
    {"Tout", "degK"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init si[] =
  {
    {"m"    , 1.},
    {"sec"  , 1.},
    {"kg"   , 1.}, 
    {"degK" , 1.},
    {"rad"  , 1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is cgs units: cm, sec, g */
struct svar_init cgs_label[] =
  {
    {"tout", "second"},
    {"lout", "cm"},
    {"aout", "cm/sec^2"},
    {"mout", "gram"},
    {"fout", "dyne"},
    {"vout", "cm/sec"},
    {"Vout", "cm^3"},
    {"dout", "g/cc"},
    {"eout", "erg"},
    {"Pout", "erg/sec"},
    {"pout", "dyne/cm^2"},
    {"Tout", "degK"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init cgs[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.},
    {"kg"   , 1000.}, 
    {"degK" ,    1.},
    {"rad"  ,    1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is cgs-ev units: cm, sec, g, eV */
struct svar_init cgs_ev_label[] =
  {
    {"tout", "second"},
    {"lout", "cm"},
    {"aout", "cm/sec^2"},
    {"mout", "gram"},
    {"fout", "dyne"},
    {"vout", "cm/sec"},
    {"Vout", "cm^3"},
    {"dout", "g/cc"},
    {"eout", "erg"},
    {"Pout", "erg/sec"},
    {"pout", "dyne/cm^2"},
    {"Tout", "eV"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init cgs_ev[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.},
    {"kg"   , 1000.}, 
    {"degK" ,    1./ 11605.},
    {"rad"  ,    1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is the shock units file: cm, usec, g */
struct svar_init shock_label[] =
  {
    {"tout", "microsecond"},
    {"lout", "cm"},
    {"aout", "cm/usec^2"},
    {"mout", "gram"},
    {"fout", "g-cm/usec^2"},
    {"vout", "cm/usec"},
    {"Vout", "cm^3"},
    {"dout", "g/cc"},
    {"eout", "g-cm^2/usec^2"},
    {"Pout", "g-cm^2/usec^3"},
    {"pout", "Mbar"},
    {"Tout", "degK"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init shock[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.0e6},
    {"kg"   , 1000.}, 
    {"degK" ,    1.},
    {"rad"  ,    1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is the "swap" units file: mm, usec, 1e-4g  */
struct svar_init swap_label[] =
  {
    {"tout", "microsecond"},
    {"lout", "mm"},
    {"aout", "mm/usec^2"},
    {"mout", "(1e-4 gram)"},
    {"fout", "(1e7 dyne)"},
    {"vout", "mm/usec"},
    {"Vout", "mm^3"},
    {"dout", "(1e-1 g/cc)"},
    {"eout", "Mega-erg"},
    {"Pout", "Mega-erg/usec"},
    {"pout", "kbar"},
    {"Tout", "degK"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init swap[] =
  {
    {"m"    ,     1000.},
    {"sec"  ,  1000000.},
    {"kg"   , 10000000.}, 
    {"degK" ,        1.},
    {"rad"  ,        1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is the ft-lbf-s units file */
struct svar_init ft_lbf_s_label[] =
  {
    {"tout", "second"},
    {"lout", "foot"},
    {"aout", "ft/sec^2"},
    {"mout", "slug"},
    {"fout", "lbf"},
    {"vout", "ft/sec"},
    {"Vout", "ft^3"},
    {"dout", "slug/ft^3"},
    {"eout", "ft-lbf"},
    {"Pout", "ft-lbf/sec"},
    {"pout", "lbf/ft^2"},
    {"Tout", "degR"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init ft_lbf_s[] =
  {
    {"m"    , 1./.3048},
    {"sec"  , 1.},
    {"kg"   , 1/4.5359237e-1/(9.806650/.3048)}, 
    {"degK" , 1.8},
    {"rad"  , 1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };


/*-------------------------------------------------------------------------------------*/
/* This is the ft-lbm-s units file */
struct svar_init ft_lbm_s_label[] =
  {
    {"tout", "second"},
    {"lout", "foot"},
    {"aout", "ft/sec^2"},
    {"mout", "lbm"},
    {"fout", "poundal"},
    {"vout", "ft/sec"},
    {"Vout", "ft^3"},
    {"dout", "lbm/ft^3"},
    {"eout", "ft-poundal"},
    {"Pout", "ft-poundal/sec"},
    {"pout", "poundal/ft^2"},
    {"Tout", "degR"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init ft_lbm_s[] =
  {
    {"m"    , 1./.3048},
    {"sec"  , 1.},
    {"kg"   , 1/.45359237}, 
    {"degK" , 1.8},
    {"rad"  , 1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

/*-------------------------------------------------------------------------------------*/
/* This is the in-lbf-s units file: inch, sec, lbf */
struct svar_init in_lbf_s_label[] =
  {
    {"tout", "second"},
    {"lout", "inch"},
    {"aout", "in/sec^2"},
    {"mout", "lbf-sec^2/in"},
    {"fout", "lbf"},
    {"vout", "in/sec"},
    {"Vout", "in^3"},
    {"dout", "lbf-sec^2/in^4"},
    {"eout", "inch-lbf"},
    {"Pout", "inch-lbf/sec"},
    {"pout", "psi"},
    {"Tout", "degR"},
    {"Aout", "radian"},
    {nullptr, nullptr}				/* Last line must be 0, 0 */
  };

struct var_init in_lbf_s[] =
  {
    {"m"    , 1./2.54e-2},
    {"sec"  , 1.},
    {"kg"   , 1/4.5359237e-1/(9.806650/2.54e-2)}, 
    {"degK" , 1.8},
    {"rad"  , 1.}, 
    {nullptr, 0}				/* Last line must be 0, 0 */
  };

struct unit_systems
{
  const char       *name;
  struct var_init  *base;
  struct svar_init *label;
};

struct unit_systems systems[] =
  {
    {"si",       si,       si_label},
    {"cgs",      cgs,      cgs_label},
    {"cgs-ev",   cgs_ev,   cgs_ev_label},
    {"shock",    shock,    shock_label},
    {"swap",     swap,     swap_label},
    {"ft-lbf-s", ft_lbf_s, ft_lbf_s_label},
    {"ft-lbm-s", ft_lbm_s, ft_lbm_s_label},
    {"in-lbf-s", in_lbf_s, in_lbf_s_label},
    {nullptr, nullptr, nullptr}
  };
} // namespace

const char *do_Units(char *type)
{
  int i, j;
  SEAMS::symrec *ptr;
  SEAMS::conv_string(type);
  
  for (i = 0; systems[i].name != nullptr; i++) {
    if (std::strcmp(type, systems[i].name) == 0) {
      break;
    }
  }
  if (systems[i].name != nullptr) {
    /* Found a match */
    for (j = 0; systems[i].label[j].vname != nullptr; j++) {
      if ((ptr = aprepro->getsym(systems[i].label[j].vname)) == nullptr) {
	ptr = aprepro->putsym(systems[i].label[j].vname, SEAMS::Aprepro::STRING_VARIABLE, true);
      }
      ptr->value.svar = systems[i].label[j].value;
    }

    for (j = 0; systems[i].base[j].vname != nullptr; j++) {
      if ((ptr = aprepro->getsym(systems[i].base[j].vname)) == nullptr) {
	ptr = aprepro->putsym(systems[i].base[j].vname, SEAMS::Aprepro::VARIABLE, true);
      }
      ptr->value.var = systems[i].base[j].value;
    }

    load_conversion(systems[i].base, systems[i].label);
    return (" ");
  }
  
    return ("Aprepro: ERROR: Invalid units system type. Valid types are: 'si', 'cgs', 'cgs-ev', 'shock', 'swap', 'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'");
  
}

void load_conversion(struct var_init *base, struct svar_init *label)
{
  SEAMS::symrec *ptr;

  const char *tout = label[ 0].value;
  const char *lout = label[ 1].value;
  const char *aout = label[ 2].value;
  const char *mout = label[ 3].value;
  const char *fout = label[ 4].value;
  const char *vout = label[ 5].value;
  const char *Vout = label[ 6].value;
  const char *dout = label[ 7].value;
  const char *eout = label[ 8].value;
  const char *Pout = label[ 9].value;
  const char *pout = label[10].value;
  const char *Tout = label[11].value;
  const char *Aout = label[12].value;

  double m    = base[0].value;
  double sec  = base[1].value;
  double kg   = base[2].value;
  double degK = base[3].value;
  double rad  = base[4].value;
  
  double foot = m * 0.3048;
  double inch = foot / 12.0;

  const char* comment = aprepro->getsym("_C_")->value.svar;
  std::string title_prefix = "\n";
  for(size_t i = 0; i < 3; i++) {
    title_prefix += comment;
}
  title_prefix += " ";

  if (echo != 0) {
    *(aprepro->infoStream)
        << title_prefix << "Outputs\n" <<
           comment << " " << std::setw(10) << "Time"   << ":\t" << tout << "\n" <<
           comment << " " << std::setw(10) << "Length" << ":\t" << lout << "\n" <<
           comment << " " << std::setw(10) << "Accel"  << ":\t" << aout << "\n" <<
           comment << " " << std::setw(10) << "Mass"   << ":\t" << mout << "\n" <<
           comment << " " << std::setw(10) << "Force"  << ":\t" << fout << "\n" <<
           comment << " " << std::setw(10) << "Velocity"  << ":\t" << vout << "\n" <<
           comment << " " << std::setw(10) << "Volume"  << ":\t" << Vout << "\n" <<
           comment << " " << std::setw(10) << "Density"  << ":\t" << dout << "\n" <<
           comment << " " << std::setw(10) << "Energy"  << ":\t" << eout << "\n" <<
           comment << " " << std::setw(10) << "Power"  << ":\t" << Pout << "\n" <<
           comment << " " << std::setw(10) << "Pressure"  << ":\t" << pout << "\n" <<
           comment << " " << std::setw(10) << "Temp"  << ":\t" << Tout << "\n" <<
           comment << " " << std::setw(10) << "Angular"  << ":\t" << Aout << "\n" <<
           comment << '\n';
  }

  if (echo != 0) {
    *(aprepro->infoStream)
        << title_prefix << "Base Dimensions\n" <<
           comment << " 1 " << std::left << std::setw(10) << "meter" << "\t= "
                   << std::setw(14) << std::setprecision(7) << m << "  " << lout << "\n" <<

           comment << " 1 " << std::left << std::setw(10) << "second" << "\t= "
                   << std::setw(14) << std::setprecision(7) << sec << "  " << tout << "\n" <<

           comment << " 1 " << std::left << std::setw(10) << "kg" << "\t= "
                   << std::setw(14) << std::setprecision(7) << kg << "  " << mout << "\n" <<

           comment << " 1 " << std::left << std::setw(10) << "kelvin" << "\t= "
                   << std::setw(14) << std::setprecision(7) << degK << "  " << Tout << "\n" <<

           comment << " 1 " << std::left << std::setw(10) << "radian" << "\t= "
                   << std::setw(14) << std::setprecision(7) << rad << "  " << Aout << '\n';
  }

  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Binary Prefixes" << '\n';
  
}DEFINE_VAR("byte", 1, "byte");
  DEFINE_VAR("KiB",  1024.0, "byte");
  DEFINE_VAR("MiB",  1024.0*1024.0, "byte");
  DEFINE_VAR("GiB",  1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("TiB",  1024.0*1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("PiB",  1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("EiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("ZiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("YiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  DEFINE_VAR("KB",   1.0e03, "byte");
  DEFINE_VAR("MB",   1.0e06, "byte");
  DEFINE_VAR("GB",   1.0e09, "byte");
  DEFINE_VAR("TB",   1.0e12, "byte");
  DEFINE_VAR("PB",   1.0e15, "byte");
  DEFINE_VAR("EB",   1.0e18, "byte");
  DEFINE_VAR("ZB",   1.0e21, "byte");
  DEFINE_VAR("YB",   1.0e24, "byte");

  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Time (T)" << '\n';
  
}DEFINE_VAR("second",                                        sec,         tout);
  DEFINE_VAR("usec",                                          sec / 1.0e6, tout);
  DEFINE_VAR("microsecond",                                   sec / 1.0e6, tout);
  DEFINE_VAR("msec",                                          sec / 1.0e3, tout);
  DEFINE_VAR("millisecond",                                   sec / 1.0e3, tout);
  DEFINE_VAR("minute",                                  60. * sec,         tout);
  DEFINE_VAR("hr",                                60. * 60. * sec,         tout);
  DEFINE_VAR("hour",                              60. * 60. * sec,         tout);
  DEFINE_VAR("day",                         24. * 60. * 60. * sec,         tout);
  DEFINE_VAR("yr",                 365.25 * 24. * 60. * 60. * sec,         tout);
  DEFINE_VAR("year",               365.25 * 24. * 60. * 60. * sec,         tout);
  DEFINE_VAR("decade",       10. * 365.25 * 24. * 60. * 60. * sec,         tout);
  DEFINE_VAR("century",     100. * 365.25 * 24. * 60. * 60. * sec,         tout);
      
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Length (L)" << '\n';
  
}DEFINE_VAR("meter",       m,         lout);
  DEFINE_VAR("metre",       m,         lout);
  DEFINE_VAR("cm",          m / 100.,  lout);
  DEFINE_VAR("centimeter",  m / 100.,  lout);
  DEFINE_VAR("centimetre",  m / 100.,  lout);
  DEFINE_VAR("mm",          m / 1000., lout);
  DEFINE_VAR("millimeter",  m / 1000., lout);
  DEFINE_VAR("millimetre",  m / 1000., lout);
  DEFINE_VAR("um",          m / 1.0e6, lout);
  DEFINE_VAR("micrometer",  m / 1.0e6, lout);
  DEFINE_VAR("micrometre",  m / 1.0e6, lout);
  DEFINE_VAR("km",          m * 1000., lout);
  DEFINE_VAR("kilometer",   m * 1000., lout);
  DEFINE_VAR("kilometre",   m * 1000., lout);
  DEFINE_VAR("ft",          foot,      lout);
  DEFINE_VAR("foot",        foot,      lout);
  DEFINE_VAR("mi",          foot * 5280., lout);
  DEFINE_VAR("mile",        foot * 5280., lout);
  DEFINE_VAR("yd",          foot * 3,  lout);
  DEFINE_VAR("yard",        foot * 3., lout);
  DEFINE_VAR("in",          inch,      lout);
  DEFINE_VAR("inch",        inch,      lout);
  DEFINE_VAR("mil",         inch / 1000., lout);

  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Acceleration (L/T^2)" << '\n';
  
}DEFINE_VAR("ga", 9.806650 * m / (sec*sec), aout);

  /* Force  (ML/T^2) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Force (ML/T^2)" << '\n';
  
}DEFINE_VAR("newton",              1.0    * kg*m/(sec*sec), fout);
  DEFINE_VAR("N",                   1.0    * kg*m/(sec*sec), fout);
  DEFINE_VAR("dyne",                1.0e-5 * kg*m/(sec*sec), fout);
  DEFINE_VAR("lbf",        4.4482216152605 * kg*m/(sec*sec), fout);
  DEFINE_VAR("kip", 1000.* 4.4482216152605 * kg*m/(sec*sec), fout);
  DEFINE_VAR("kgf",               9.806650 * kg*m/(sec*sec), fout);
  DEFINE_VAR("gf",                9.806650 * kg*m/(sec*sec)/1000., fout);
  DEFINE_VAR("pdl",            1.382550e-1 * kg*m/(sec*sec), fout);
  DEFINE_VAR("poundal",        1.382550e-1 * kg*m/(sec*sec), fout);
  DEFINE_VAR("ounce",      4.4482216152605 * kg*m/(sec*sec)/16.0, fout);

  /* Mass (M) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Mass (M)" << '\n';
  
}DEFINE_VAR("gram",              kg / 1000., mout);
  DEFINE_VAR("g",                 kg / 1000., mout);
  DEFINE_VAR("lbm",   453.59237 * kg / 1000., mout);
  DEFINE_VAR("slug",  453.59237 * kg / 1000. * 32.17404856, mout);
  DEFINE_VAR("lbfs2pin",  4.4482216152605 * kg/0.0254, mout);
  
  /* Velocity (L/T) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Velocity (L/T)" << '\n';
  
}DEFINE_VAR("mps",  m/sec, vout);
  DEFINE_VAR("fps",  foot / sec, vout);
  DEFINE_VAR("mph",  (foot * 5280.) / (60. * 60. * sec), vout);
  DEFINE_VAR("ips",  (inch) / sec, vout);
  DEFINE_VAR("kph",  (1000. * m) / (60. * 60. * sec), vout);
  DEFINE_VAR("kps",  (1000. * m) / sec, vout);

  /* Volume (L^3) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Volume (L^3)" << '\n';
  
}DEFINE_VAR("liter",              (m*m*m)/1000., Vout);
  DEFINE_VAR("gal",     3.785412 * (m*m*m)/1000., Vout);
  DEFINE_VAR("gallon",  3.785412 * (m*m*m)/1000., Vout);

  /* Density (M/L^3) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Density (M/L^3)" << '\n';
  
}DEFINE_VAR("gpcc",      (kg/1000.)/((m/100.)*(m/100.)*(m/100.)), dout);
  DEFINE_VAR("kgpm3",      kg /(m*m*m), dout);
  DEFINE_VAR("lbfs2pin4",  (4.4482216152605 * kg*m/(sec*sec))*sec*sec / (inch*inch*inch*inch), dout);
  DEFINE_VAR("lbmpin3",    (453.59237 * kg / 1000.) / (inch*inch*inch), dout);
  DEFINE_VAR("lbmpft3",    (453.59237 * kg / 1000.) / (foot*foot*foot), dout);
  DEFINE_VAR("slugpft3",   (453.59237 * kg / 1000. * 32.17404856) / (foot*foot*foot), dout);

  /* Power: (M L^2 / T^3) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Power (M L^2 / T^3)" << '\n';
  
}DEFINE_VAR("W",    kg*m/(sec*sec)*m/sec, Pout);
  DEFINE_VAR("watt", kg*m/(sec*sec)*m/sec, Pout);
  DEFINE_VAR("Hp",   kg*m/(sec*sec)*m/sec * 746, Pout); /* --- (electric horsepower) */

  /* Energy (ML^2/T^2) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Energy (M L^2 / T^2)" << '\n';
  
}DEFINE_VAR("joule",  kg*m/(sec*sec)*m, eout);
  DEFINE_VAR("J",      kg*m/(sec*sec)*m, eout);
  DEFINE_VAR("ftlbf",  kg*m/(sec*sec)*m * 1.355818, eout); 
  DEFINE_VAR("Btu",    kg*m/(sec*sec)*m * 1.05505585262e3, eout); /*--- I18n Table */
  DEFINE_VAR("erg",    kg*m/(sec*sec)*m * 1.0e-7, eout);
  DEFINE_VAR("calorie",kg*m/(sec*sec)*m * 4.18680, eout);  /* --- I18n Table */
  DEFINE_VAR("kwh",    kg*m/(sec*sec)*m * 1000.0 * 60. * 60., eout); 
  DEFINE_VAR("therm",  kg*m/(sec*sec)*m * 1.054804e8, eout); /*       --- U.S. */
  DEFINE_VAR("tonTNT",  kg*m/(sec*sec)*m * 4.184e9, eout); 

  /* Pressure: (M/L/T^2) */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Pressure (M/L/T^2)" << '\n';
  
}DEFINE_VAR("Pa",      kg*m/(sec*sec) / (m*m), pout);
  DEFINE_VAR("pascal",  kg*m/(sec*sec) / (m*m), pout);
  DEFINE_VAR("MPa",     kg*m/(sec*sec) / (m*m) * 1.0e6, pout); 
  DEFINE_VAR("GPa",     kg*m/(sec*sec) / (m*m) * 1.0e9, pout); 
  DEFINE_VAR("bar",     kg*m/(sec*sec) / (m*m) * 1.0e5, pout); 
  DEFINE_VAR("kbar",    kg*m/(sec*sec) / (m*m) * 1.0e5 * 1.0e3, pout); 
  DEFINE_VAR("Mbar",    kg*m/(sec*sec) / (m*m) * 1.0e5 * 1.0e6, pout); 
  DEFINE_VAR("psi",     4.4482216152605 * kg*m/(sec*sec) / (inch*inch), pout); 
  DEFINE_VAR("ksi",     4.4482216152605 * kg*m/(sec*sec) / (inch*inch) * 1000.0, pout); 
  DEFINE_VAR("psf",     4.4482216152605 * kg*m/(sec*sec) / (foot*foot), pout); 
  DEFINE_VAR("atm",     kg*m/(sec*sec) / (m*m) * 1.013250e5, pout);  /* --- std atmosphere */
  DEFINE_VAR("torr",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0, pout); 
  DEFINE_VAR("mHg",     kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0 * 1000.0, pout); 
  DEFINE_VAR("mmHg",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0, pout); 
  DEFINE_VAR("inHg",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0 * 25.4, pout);
  DEFINE_VAR("inH2O",   kg*m/(sec*sec) / (m*m) * 249.082, pout); 
  DEFINE_VAR("ftH2O",   kg*m/(sec*sec) / (m*m) * 249.082 * 12.0, pout); 

  /* Temperature: */
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Temperature" << '\n';
  
}DEFINE_VAR("kelvin",         degK, Tout); 
  DEFINE_VAR("degC",           degK, Tout); 
  DEFINE_VAR("degF",   5./9. * degK, Tout); 
  DEFINE_VAR("degR",   5./9. * degK, Tout);
  DEFINE_VAR("rankine",5./9. * degK, Tout);
  DEFINE_VAR("eV",     11605 * degK, Tout);

  /* Angular */
#define PI  3.141592653589793238462643
  if (echo != 0) { *(aprepro->infoStream) << title_prefix << "Angular" << '\n';
  
}DEFINE_VAR("rev",    2.0 * PI * rad, Aout); 
  DEFINE_VAR("deg",    2.0 * PI * rad / 360.0, Aout); 
  DEFINE_VAR("degree", 2.0 * PI * rad / 360.0, Aout); 
  DEFINE_VAR("arcmin", 2.0 * PI * rad / 360.0 / 60.0, Aout); 
  DEFINE_VAR("arcsec", 2.0 * PI * rad / 360.0 / 60.0 / 60.0, Aout); 
  DEFINE_VAR("grade",  2.0 * PI * rad / 360.0 * 0.9, Aout); 
}


}  // namespace SEAMS
