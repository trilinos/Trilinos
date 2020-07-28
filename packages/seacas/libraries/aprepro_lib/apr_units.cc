// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

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
  extern bool     echo;

  // clang-format off
namespace {
  std::string &comment()
  {
    static std::string com = aprepro->getsym("_C_")->value.svar;
    return com;
  }

  void define_var(const char *name, double val, const char *label)
  {
    aprepro->add_variable(name, val, true);
    if (echo) {
      *(aprepro->infoStream) << comment() << " 1 " << std::left << std::setw(10) << name
                             << "\t= " << std::setw(14) << std::setprecision(7) << val << "  "
                             << label << '\n';
    }
  }

  void load_conversion(var_init *base, svar_init *label);
  constexpr double LBF_TO_N = 4.4482216152605;
  constexpr double PI = 3.141592653589793238462643;

//-------------------------------------------------------------------------------------
// SI Units
svar_init si_label[] =
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
    {nullptr, nullptr}
  };

var_init si[] =
  {
    {"m"    , 1.},
    {"sec"  , 1.},
    {"kg"   , 1.},
    {"degK" , 1.},
    {"rad"  , 1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is cgs units: cm, sec, g
svar_init cgs_label[] =
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
    {nullptr, nullptr}
  };

var_init cgs[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.},
    {"kg"   , 1000.},
    {"degK" ,    1.},
    {"rad"  ,    1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is cgs-ev units: cm, sec, g, eV
svar_init cgs_ev_label[] =
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
    {nullptr, nullptr}
  };

var_init cgs_ev[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.},
    {"kg"   , 1000.},
    {"degK" ,    1./ 11604.5221},
    {"rad"  ,    1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is the shock units file: cm, usec, g
svar_init shock_label[] =
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
    {nullptr, nullptr}
  };

var_init shock[] =
  {
    {"m"    ,  100.},
    {"sec"  ,    1.0e6},
    {"kg"   , 1000.},
    {"degK" ,    1.},
    {"rad"  ,    1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is the "swap" units file: mm, usec, 1e-4g
svar_init swap_label[] =
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
    {nullptr, nullptr}
  };

var_init swap[] =
  {
    {"m"    ,     1000.},
    {"sec"  ,  1000000.},
    {"kg"   , 10000000.},
    {"degK" ,        1.},
    {"rad"  ,        1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is the ft-lbf-s units file
svar_init ft_lbf_s_label[] =
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
    {nullptr, nullptr}
  };

var_init ft_lbf_s[] =
  {
    {"m"    , 1./.3048},
    {"sec"  , 1.},
    {"kg"   , 1/4.5359237e-1/(9.806650/.3048)},
    {"degK" , 1.8},
    {"rad"  , 1.},
    {nullptr, 0}
  };


//-------------------------------------------------------------------------------------
// This is the ft-lbm-s units file
svar_init ft_lbm_s_label[] =
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
    {nullptr, nullptr}
  };

var_init ft_lbm_s[] =
  {
    {"m"    , 1./.3048},
    {"sec"  , 1.},
    {"kg"   , 1/.45359237},
    {"degK" , 1.8},
    {"rad"  , 1.},
    {nullptr, 0}
  };

//-------------------------------------------------------------------------------------
// This is the in-lbf-s units file: inch, sec, lbf
svar_init in_lbf_s_label[] =
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
    {nullptr, nullptr}
  };

var_init in_lbf_s[] =
  {
    {"m"    , 1./2.54e-2},
    {"sec"  , 1.},
    {"kg"   , 1/4.5359237e-1/(9.806650/2.54e-2)},
    {"degK" , 1.8},
    {"rad"  , 1.},
    {nullptr, 0}
  };

struct unit_systems
{
  const char       *name;
  var_init  *base;
  svar_init *label;
};

unit_systems systems[] =
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
    SEAMS::conv_string(type);

    int i;
    for (i = 0; systems[i].name != nullptr; i++) {
      if (std::strcmp(type, systems[i].name) == 0) {
        break;
      }
    }
    if (systems[i].name != nullptr) {
      // Found a match
      for (int j = 0; systems[i].label[j].vname != nullptr; j++) {
        aprepro->add_variable(systems[i].label[j].vname, systems[i].label[j].value, true);
      }

      for (int j = 0; systems[i].base[j].vname != nullptr; j++) {
        aprepro->add_variable(systems[i].base[j].vname, systems[i].base[j].value, true);
      }

      load_conversion(systems[i].base, systems[i].label);
      return (" ");
    }

    return ("Aprepro: ERROR: Invalid units system type. Valid types are: 'si', 'cgs', 'cgs-ev', "
            "'shock', 'swap', 'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'");
  }

namespace {
void load_conversion(var_init *base, svar_init *label)
{
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

  std::string title_prefix = "\n";
  for(size_t i = 0; i < 3; i++) {
    title_prefix += comment();
}
  title_prefix += " ";

  if (echo != 0) {
    *(aprepro->infoStream)
        << title_prefix << "Outputs\n" <<
           comment() << " " << std::setw(10) << "Time"   << ":\t" << tout << "\n" <<
           comment() << " " << std::setw(10) << "Length" << ":\t" << lout << "\n" <<
           comment() << " " << std::setw(10) << "Accel"  << ":\t" << aout << "\n" <<
           comment() << " " << std::setw(10) << "Mass"   << ":\t" << mout << "\n" <<
           comment() << " " << std::setw(10) << "Force"  << ":\t" << fout << "\n" <<
           comment() << " " << std::setw(10) << "Velocity"  << ":\t" << vout << "\n" <<
           comment() << " " << std::setw(10) << "Volume"  << ":\t" << Vout << "\n" <<
           comment() << " " << std::setw(10) << "Density"  << ":\t" << dout << "\n" <<
           comment() << " " << std::setw(10) << "Energy"  << ":\t" << eout << "\n" <<
           comment() << " " << std::setw(10) << "Power"  << ":\t" << Pout << "\n" <<
           comment() << " " << std::setw(10) << "Pressure"  << ":\t" << pout << "\n" <<
           comment() << " " << std::setw(10) << "Temp"  << ":\t" << Tout << "\n" <<
           comment() << " " << std::setw(10) << "Angular"  << ":\t" << Aout << "\n" <<
           comment() << '\n';
  }

  if (echo != 0) {
    *(aprepro->infoStream)
        << title_prefix << "Base Dimensions\n" <<
           comment() << " 1 " << std::left << std::setw(10) << "meter" << "\t= "
                   << std::setw(14) << std::setprecision(7) << m << "  " << lout << "\n" <<

           comment() << " 1 " << std::left << std::setw(10) << "second" << "\t= "
                   << std::setw(14) << std::setprecision(7) << sec << "  " << tout << "\n" <<

           comment() << " 1 " << std::left << std::setw(10) << "kg" << "\t= "
                   << std::setw(14) << std::setprecision(7) << kg << "  " << mout << "\n" <<

           comment() << " 1 " << std::left << std::setw(10) << "kelvin" << "\t= "
                   << std::setw(14) << std::setprecision(7) << degK << "  " << Tout << "\n" <<

           comment() << " 1 " << std::left << std::setw(10) << "radian" << "\t= "
                   << std::setw(14) << std::setprecision(7) << rad << "  " << Aout << '\n';
  }

  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Binary Prefixes" << '\n';
  }
  define_var("byte", 1, "byte");
  define_var("KiB",  1024.0, "byte");
  define_var("MiB",  1024.0*1024.0, "byte");
  define_var("GiB",  1024.0*1024.0*1024.0, "byte");
  define_var("TiB",  1024.0*1024.0*1024.0*1024.0, "byte");
  define_var("PiB",  1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  define_var("EiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  define_var("ZiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  define_var("YiB",  1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0*1024.0, "byte");
  define_var("KB",   1.0e03, "byte");
  define_var("MB",   1.0e06, "byte");
  define_var("GB",   1.0e09, "byte");
  define_var("TB",   1.0e12, "byte");
  define_var("PB",   1.0e15, "byte");
  define_var("EB",   1.0e18, "byte");
  define_var("ZB",   1.0e21, "byte");
  define_var("YB",   1.0e24, "byte");

  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Time (T)" << '\n';
  }
  define_var("second",                                        sec,         tout);
  define_var("usec",                                          sec / 1.0e6, tout);
  define_var("microsecond",                                   sec / 1.0e6, tout);
  define_var("msec",                                          sec / 1.0e3, tout);
  define_var("millisecond",                                   sec / 1.0e3, tout);
  define_var("minute",                                  60. * sec,         tout);
  define_var("hr",                                60. * 60. * sec,         tout);
  define_var("hour",                              60. * 60. * sec,         tout);
  define_var("day",                         24. * 60. * 60. * sec,         tout);
  define_var("yr",                 365.25 * 24. * 60. * 60. * sec,         tout);
  define_var("year",               365.25 * 24. * 60. * 60. * sec,         tout);
  define_var("decade",       10. * 365.25 * 24. * 60. * 60. * sec,         tout);
  define_var("century",     100. * 365.25 * 24. * 60. * 60. * sec,         tout);

  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Length (L)" << '\n';
  }
  define_var("meter",       m,         lout);
  define_var("metre",       m,         lout);
  define_var("cm",          m / 100.,  lout);
  define_var("centimeter",  m / 100.,  lout);
  define_var("centimetre",  m / 100.,  lout);
  define_var("mm",          m / 1000., lout);
  define_var("millimeter",  m / 1000., lout);
  define_var("millimetre",  m / 1000., lout);
  define_var("um",          m / 1.0e6, lout);
  define_var("micrometer",  m / 1.0e6, lout);
  define_var("micrometre",  m / 1.0e6, lout);
  define_var("km",          m * 1000., lout);
  define_var("kilometer",   m * 1000., lout);
  define_var("kilometre",   m * 1000., lout);
  define_var("ft",          foot,      lout);
  define_var("foot",        foot,      lout);
  define_var("mi",          foot * 5280., lout);
  define_var("mile",        foot * 5280., lout);
  define_var("yd",          foot * 3,  lout);
  define_var("yard",        foot * 3., lout);
  define_var("in",          inch,      lout);
  define_var("inch",        inch,      lout);
  define_var("mil",         inch / 1000., lout);

  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Acceleration (L/T^2)" << '\n';
  }
  define_var("ga", 9.806650 * m / (sec*sec), aout);

  // Force  (ML/T^2)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Force (ML/T^2)" << '\n';
  }
  define_var("newton",              1.0    * kg*m/(sec*sec), fout);
  define_var("N",                   1.0    * kg*m/(sec*sec), fout);
  define_var("dyne",                1.0e-5 * kg*m/(sec*sec), fout);
  define_var("lbf",        LBF_TO_N * kg*m/(sec*sec), fout);
  define_var("kip", 1000.* LBF_TO_N * kg*m/(sec*sec), fout);
  define_var("kgf",               9.806650 * kg*m/(sec*sec), fout);
  define_var("gf",                9.806650 * kg*m/(sec*sec)/1000., fout);
  define_var("pdl",            1.382550e-1 * kg*m/(sec*sec), fout);
  define_var("poundal",        1.382550e-1 * kg*m/(sec*sec), fout);
  define_var("ounce",      LBF_TO_N * kg*m/(sec*sec)/16.0, fout);

  // Mass (M)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Mass (M)" << '\n';
  }
  define_var("gram",              kg / 1000., mout);
  define_var("g",                 kg / 1000., mout);
  define_var("lbm",   453.59237 * kg / 1000., mout);
  define_var("slug",  453.59237 * kg / 1000. * 32.17404856, mout);
  define_var("lbfs2pin",  LBF_TO_N * kg/0.0254, mout);

  // Velocity (L/T)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Velocity (L/T)" << '\n';
  }
  define_var("mps",  m/sec, vout);
  define_var("fps",  foot / sec, vout);
  define_var("mph",  (foot * 5280.) / (60. * 60. * sec), vout);
  define_var("ips",  (inch) / sec, vout);
  define_var("kph",  (1000. * m) / (60. * 60. * sec), vout);
  define_var("kps",  (1000. * m) / sec, vout);

  // Volume (L^3)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Volume (L^3)" << '\n';
  }
  define_var("liter",              (m*m*m)/1000., Vout);
  define_var("gal",     3.785412 * (m*m*m)/1000., Vout);
  define_var("gallon",  3.785412 * (m*m*m)/1000., Vout);

  // Density (M/L^3)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Density (M/L^3)" << '\n';
  }
  define_var("gpcc",      (kg/1000.)/((m/100.)*(m/100.)*(m/100.)), dout);
  define_var("kgpm3",      kg /(m*m*m), dout);
  define_var("lbfs2pin4",  (LBF_TO_N * kg*m/(sec*sec))*sec*sec / (inch*inch*inch*inch), dout);
  define_var("lbmpin3",    (453.59237 * kg / 1000.) / (inch*inch*inch), dout);
  define_var("lbmpft3",    (453.59237 * kg / 1000.) / (foot*foot*foot), dout);
  define_var("slugpft3",   (453.59237 * kg / 1000. * 32.17404856) / (foot*foot*foot), dout);

  // Power: (M L^2 / T^3)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Power (M L^2 / T^3)" << '\n';
  }
  define_var("W",    kg*m/(sec*sec)*m/sec, Pout);
  define_var("watt", kg*m/(sec*sec)*m/sec, Pout);
  define_var("Hp",   kg*m/(sec*sec)*m/sec * 746, Pout); // --- (electric horsepower)

  // Energy (ML^2/T^2)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Energy (M L^2 / T^2)" << '\n';
  }
  define_var("joule",   kg*m/(sec*sec)*m, eout);
  define_var("J",       kg*m/(sec*sec)*m, eout);
  define_var("ftlbf",   kg*m/(sec*sec)*m * 1.355818, eout);
  define_var("Btu",     kg*m/(sec*sec)*m * 1.05505585262e3, eout); //--- I18n Table
  define_var("erg",     kg*m/(sec*sec)*m * 1.0e-7, eout);
  define_var("calorie", kg*m/(sec*sec)*m * 4.18680, eout);  // --- I18n Table
  define_var("kwh",     kg*m/(sec*sec)*m * 1000.0 * 60. * 60., eout);
  define_var("therm",   kg*m/(sec*sec)*m * 1.054804e8, eout); //       --- U.S.
  define_var("tonTNT",  kg*m/(sec*sec)*m * 4.184e9, eout);

  // Pressure: (M/L/T^2)
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Pressure (M/L/T^2)" << '\n';
  }
  define_var("Pa",      kg*m/(sec*sec) / (m*m), pout);
  define_var("pascal",  kg*m/(sec*sec) / (m*m), pout);
  define_var("MPa",     kg*m/(sec*sec) / (m*m) * 1.0e6, pout);
  define_var("GPa",     kg*m/(sec*sec) / (m*m) * 1.0e9, pout);
  define_var("bar",     kg*m/(sec*sec) / (m*m) * 1.0e5, pout);
  define_var("kbar",    kg*m/(sec*sec) / (m*m) * 1.0e5 * 1.0e3, pout);
  define_var("Mbar",    kg*m/(sec*sec) / (m*m) * 1.0e5 * 1.0e6, pout);
  define_var("psi",     LBF_TO_N * kg*m/(sec*sec) / (inch*inch), pout);
  define_var("ksi",     LBF_TO_N * kg*m/(sec*sec) / (inch*inch) * 1000.0, pout);
  define_var("psf",     LBF_TO_N * kg*m/(sec*sec) / (foot*foot), pout);
  define_var("atm",     kg*m/(sec*sec) / (m*m) * 1.013250e5, pout);  // --- std atmosphere
  define_var("torr",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0, pout);
  define_var("mHg",     kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0 * 1000.0, pout);
  define_var("mmHg",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0, pout);
  define_var("inHg",    kg*m/(sec*sec) / (m*m) * 1.013250e5 / 760.0 * 25.4, pout);
  define_var("inH2O",   kg*m/(sec*sec) / (m*m) * 249.082, pout);
  define_var("ftH2O",   kg*m/(sec*sec) / (m*m) * 249.082 * 12.0, pout);

  // Temperature:
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Temperature" << '\n';
  }
  define_var("kelvin",         degK, Tout);
  define_var("degC",           degK, Tout);
  define_var("degF",   5./9. * degK, Tout);
  define_var("degR",   5./9. * degK, Tout);
  define_var("rankine",5./9. * degK, Tout);
  define_var("eV",     11604.5221 * degK, Tout);

  // Angular
  if (echo != 0) {
    *(aprepro->infoStream) << title_prefix << "Angular" << '\n';
  }
  define_var("rev",    2.0 * PI * rad, Aout);
  define_var("deg",    2.0 * PI * rad / 360.0, Aout);
  define_var("degree", 2.0 * PI * rad / 360.0, Aout);
  define_var("arcmin", 2.0 * PI * rad / 360.0 / 60.0, Aout);
  define_var("arcsec", 2.0 * PI * rad / 360.0 / 60.0 / 60.0, Aout);
  define_var("grade",  2.0 * PI * rad / 360.0 * 0.9, Aout);
}
}
}  // namespace SEAMS
