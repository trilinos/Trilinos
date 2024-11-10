// Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
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
    aprepro->add_variable(name, val, true, true);
    if (echo) {
      *(aprepro->infoStream) << comment() << " 1 " << std::left << std::setw(10) << name
                             << "\t= " << std::setw(14) << std::setprecision(7) << val << "  "
                             << label << '\n';
    }
  }

  void load_conversion(var_init *base, svar_init *label);
  void load_constants(var_init *base);
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
    {"m"    ,      1'000.},
    {"sec"  ,  1'000'000.},
    {"kg"   , 10'000'000.},
    {"degK" ,          1.},
    {"rad"  ,          1.},
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
        aprepro->add_variable(systems[i].label[j].vname, systems[i].label[j].value, true, true);
      }

      for (int j = 0; systems[i].base[j].vname != nullptr; j++) {
        aprepro->add_variable(systems[i].base[j].vname, systems[i].base[j].value, true, true);
      }

      symrec *var = aprepro->getsym("_UNITS_SYSTEM");
      var->value.svar = type;

      load_conversion(systems[i].base, systems[i].label);
      load_constants(systems[i].base);
      return (" ");
    }

    return ("Aprepro: ERROR: Invalid units system type. Valid types are: 'si', 'cgs', 'cgs-ev', "
            "'shock', 'swap', 'ft-lbf-s', 'ft-lbm-s', 'in-lbf-s'");
  }

namespace {
  void load_constants(var_init *base)
  {
    double Ohm     = 1;
    double Hz      = 1;
    double A  = 1;
    double C  = 1;
    double F   = 1;
    double V    = 1;
    double T   = 1;
    double mol    = 1;
    double S = 1;
    double GeV     = 1;

    double m      = base[0].value;
    double sec    = base[1].value;
    double kg     = base[2].value;
    double degK   = base[3].value;

    double N = kg * m / (sec * sec);
    double J = N * m;
    double W = J / sec;

    double sec2   = sec * sec;
    double degK4  = degK * degK * degK * degK;
    double m2     = m * m;
    double m3     = m * m * m;
    double C2     = C * C;
    double A2     = A * A;

    if (echo != 0) {
      std::string title_prefix = "\n";
      for (size_t i = 0; i < 3; i++) {
	title_prefix += comment();
      }
      title_prefix += " ";
      *(aprepro->infoStream)
	<< title_prefix
	<< "Physical Constants (https://en.wikipedia.org/wiki/List_of_physical_constants)"
	<< '\n';
    }
    define_var("Avogadro_constant",                 6.02214076E23 / mol,               "# 6.02214076E23 / mol");
    define_var("Bohr_magneton",                     9.2740100783E-24 * J / T,          "# 9.2740100783E-24 J/T");
    define_var("Bohr_radius",                       5.29177210903E-11 * m,             "# 5.29177210903E-11 m");
    define_var("Boltzmann_constant",                1.380649E-23 * J / degK,           "# 1.380649E-23 J/degK");
    define_var("Coulomb_constant",                  8.9875517923E9 * N * m2 / C2,      "# 8.9875517923E9 N m^2 C^-2");
    define_var("Faraday_constant",                  96485.3321233100184 * C / mol,     "# 96485.3321233100184 C/mol");
    define_var("Fermi_coupling_constant",           1.166378710E-5 / (GeV * GeV),      "# 1.166378710E-5 GeV^-2");
    define_var("Hartree_energy",                    4.3597447222071E-18 * J,           "# 4.3597447222071E-18 J");
    define_var("Josephson_constant",                483597.8484E9 * Hz / V,            "# 483597.8484E9 Hz/V");
    define_var("Newtonian_constant_of_gravitation", 6.67430E-11 * m3 / kg / sec2,      "# 6.67430E-11 m^3/kg s^-2");
    define_var("Gravitational_constant",            6.67430E-11 * m3 / kg / sec2,      "# 6.67430E-11 m^3/kg s^-2");
    define_var("Planck_constant",                   6.62607015E-34 * J / Hz,           "# 6.62607015E-34 J/Hz");
    define_var("Rydberg_constant",                  10973731.568160 / m,               "# 10973731.568160 m");
    define_var("Rydberg_unit_of_energy",            2.1798723611035E-18 * J,           "# 2.1798723611035E-18 J");
    define_var("Stefan_Boltzmann_constant",         5.670374419E-8 * W / m2 / degK4,   "# 5.670374419E-8 W m^-2 degK^-4");
    define_var("Thomson_cross_section",             6.6524587321E-29 * m2,             "# 6.6524587321E-29 m^2");
    define_var("W_to_Z_mass_ratio",                 0.88153,                           "");
    define_var("Wien_entropy_displacement_law_constant",    3.002916077E-3 * m * degK, "# 3.002916077E-3 m degK");
    define_var("Wien_frequency_displacement_law_constant",  5.878925757E10 * Hz / degK,"# 5.878925757E10 Hz/degK");
    define_var("Wien_wavelength_displacement_law_constant", 2.897771955E-3 * m * degK, "# 2.897771955E-3 m degK");
    define_var("atomic_mass_constant",              1.66053906660E-27 * kg,            "# 1.66053906660E-27 kg");
    define_var("atomic_mass_of_carbon_12",          1.99264687992E-26 * kg,            "# 1.99264687992E-26 kg");
    define_var("characteristic_impedance_of_vacuum", 376.730313668 * Ohm,              "# 376.730313668 Ohm");
    define_var("classical_electron_radius",         2.8179403262E-15 * m,              "# 2.8179403262E-15 m");
    define_var("conductance_quantum",               7.748091729E-5 * S,                "# 7.748091729E-5 S");
    define_var("cosmological_constant",             1.089E-52 / m2,                    "# 1.089E-52 m^-2");
    define_var("electron_g_factor",                 -2.00231930436256,                 "");
    define_var("electron_mass",                     9.1093837015E-31 * kg,             "# 9.1093837015E-31 kg");
    define_var("elementary_charge",                 1.602176634E-19 * C,               "# 1.602176634E-19 C");
    define_var("fine_structure_constant",           7.2973525693E-3,                   "");
    define_var("first_radiation_constant",          3.741771852E-16 * W * m2,          "# 3.741771852E-16 W m^2");
    define_var("hyperfine_transition_frequency_of_133Cs",   9192631770 * Hz,           "# 9192631770 Hz");
    define_var("inverse_conductance_quantum",       12906.40372 * Ohm,                 "# 12906.40372 Ohm");
    define_var("inverse_fine_structure_constant",   137.035999084,                     "");
    define_var("magnetic_flux_quantum",             2.067833848E-15 * V * sec,         "# 2.067833848E-15 V s");
    define_var("molar_Planck_constant",             3.9903127128934314E-10 * J * sec / mol, "# 3.9903127128934314E-10 J s / mol");
    define_var("molar_gas_constant",                8.31446261815324 * J / mol / degK, "# 8.31446261815324 J/mol/degK");
    define_var("molar_mass_constant",               0.99999999965E-3 * kg / mol,       "# 0.99999999965E-3 kg/mol");
    define_var("molar_mass_of_carbon_12",           11.9999999958E-3 * kg / mol,       "# 11.9999999958E-3 kg/mol");
    define_var("muon_g_factor",                     -2.0023318418,                     "");
    define_var("muon_mass",                         1.883531627E-28 * kg,              "# 1.883531627E-28 kg");
    define_var("neutron_mass",                      1.67492749804E-27 * kg,            "# 1.67492749804E-27 kg");
    define_var("nuclear_magneton",                  5.0507837461E-27 * J / T,          "# 5.0507837461E-27  J/T");
    define_var("proton_g_factor",                   5.5856946893,                      "");
    define_var("proton_mass",                       1.67262192369E-27 * kg,            "# 1.67262192369E-27 kg");
    define_var("proton_to_electron_mass_ratio",     1836.15267343,                     "");
    define_var("quantum_of_circulation",            3.6369475516E-4 * m2 / sec,        "# 3.6369475516E-4 m^2/s");
    define_var("reduced_Planck_constant",           1.054571817E-34 * J * sec,         "# 1.054571817E-34 J s");
    define_var("sec_radiation_constant",            1.438776877E-2 * m * degK,         "# 1.438776877E-2 m degK");
    define_var("speed_of_light_in_vacuum",          299792458 * m / sec,               "# 299792458 m/s");
    define_var("tau_mass",                          3.16754E-27 * kg,                  "# 3.16754E-27 kg");
    define_var("top_quark_mass",                    3.0784E-25 * kg,                   "# 3.0784E-25 kg");
    define_var("vacuum_electric_permittivity",      8.8541878128E-12 * F / m,          "# 8.8541878128E-12 F/m");
    define_var("vacuum_magnetic_permeability",      1.25663706212E-6 * N / A2,         "# 1.25663706212E-6 N A^-2");
    define_var("von_Klitzing_constant",             25812.80745 * Ohm,                 "# 25812.80745 Ohm");
    define_var("weak_mixing_angle",                 0.22290,                           "");
  }

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
}  // namespace
}  // namespace SEAMS
