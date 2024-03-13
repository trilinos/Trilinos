// Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "aprepro.h" // for Aprepro
#include <iostream>  // for cout, ostream, etc
#include <string>    // for string, operator<<
#include <vector>    // for vector

// This function is used below in the example showing how an
// application can add its own functions to an aprepro instance.
double succ(double i) { return ++i; }

int main(int argc, char *argv[])
{
  SEAMS::Aprepro aprepro;

  // Read and parse a string's worth of data at a time.
  // Cannot use looping/ifs/... with this method.
  std::vector<std::string> strings;
  strings.emplace_back("$ Test program for Aprepro");
  strings.emplace_back("$");
  strings.emplace_back("Test number representations");
  strings.emplace_back("{1}        {10e-1} {10.e-1}        {.1e+1} {.1e1}");
  strings.emplace_back("{1}        {10E-1} {10.E-1}        {.1E+1} {.1E1}");
  strings.emplace_back("   ");
  strings.emplace_back("Test assign statements:");
  strings.emplace_back("{_a = 5}   {b=_a}  $ Should print 5 5");
  strings.emplace_back("{_a +=b}   {_a}    $ Should print 10 10");
  strings.emplace_back("{_a -=b}   {_a}    $ Should print 5 5");
  strings.emplace_back("{_a *=b}   {_a}    $ Should print 25 25");
  strings.emplace_back("{_a /=b}   {_a}    $ Should print 5 5");
  strings.emplace_back("{_a ^=b}   {_a}    $ Should print 3125 3125");
  strings.emplace_back("{_a = b}   {_a**=b}        {_a}    $ Should print 5 3125 3125");
  strings.emplace_back("");
  strings.emplace_back("Test trigonometric functions (radians)");
  strings.emplace_back("{pi = d2r(180)} {atan2(0,-1)} {4*atan(1.0)} $ Three values of pi");
  strings.emplace_back("{_a = sin(pi/4)}   {pi-4*asin(_a)} $ sin(pi/4)");
  strings.emplace_back("{_b = cos(pi/4)}   {pi-4*acos(_b)} $ cos(pi/4)");
  strings.emplace_back("{_c = tan(pi/4)}   {pi-4*atan(_c)} $ tan(pi/4)");
  strings.emplace_back("");
  strings.emplace_back("Test trigonometric functions (degrees)");
  strings.emplace_back("{r2d(pi)}  {pid = atan2d(0,-1)}    {4 * atand(1.0)}");
  strings.emplace_back("{ad = sind(180/4)} {180-4*asind(ad)}       $ sin(180/4)");
  strings.emplace_back("{bd = cosd(180/4)} {180-4*acosd(bd)}       $ cos(180/4)");
  strings.emplace_back("{cd = tand(180/4)} {180-4*atand(cd)}       $ tan(180/4)");
  strings.emplace_back("");
  strings.emplace_back("Test max, min, sign, dim, abs");
  strings.emplace_back("{pmin = min(0.5, 1.0)}     {nmin = min(-0.5, -1.0)} $ Should be 0.5, -1");
  strings.emplace_back("{pmax = max(0.5, 1.0)}     {nmax = max(-0.5, -1.0)} $ Should be 1.0, -0.5");
  strings.emplace_back("{zero = 0} {sign(0.5, zero) + sign(0.5, -zero)}    $ Should be 0 1");
  strings.emplace_back("{nonzero = 1} {sign(0.5, nonzero) + sign(0.5, -nonzero)} $ Should be 1 0");
  strings.emplace_back("{dim(5.5, 4.5)}    {dim(4.5, 5.5)} $ Should be 1 0");
  strings.emplace_back("");
  strings.emplace_back("{ifyes = 1} {ifno = 0}");
  strings.emplace_back("$ Test ternary...");
  strings.emplace_back(R("{ifyes == 1 ? " correct " : " incorrect "}"));
  strings.emplace_back(R("{ifno == 1 ? " incorrect " : " correct "}"));
  strings.emplace_back("");
  strings.emplace_back("$ Test ifdef lines");
  strings.emplace_back("   {Ifdef(ifyes)}");
  strings.emplace_back("This line should be echoed. (a)");
  strings.emplace_back(" {Endif}");
  strings.emplace_back("This line should be echoed. (b)");
  strings.emplace_back("     {Ifdef(ifno)}");
  strings.emplace_back("This line should not be echoed");
  strings.emplace_back("    {Endif}");
  strings.emplace_back("This line should be echoed. (c)");
  strings.emplace_back("  {Ifdef(ifundefined)}");
  strings.emplace_back("This line should not be echoed");
  strings.emplace_back("        {Endif}");
  strings.emplace_back("This line should be echoed. (d)");
  strings.emplace_back("");
  strings.emplace_back("$ Test ifdef - else lines");
  strings.emplace_back("             {Ifdef(ifyes)}");
  strings.emplace_back("This line should be echoed. (1)");
  strings.emplace_back("                   {Else}");
  strings.emplace_back("This line should not be echoed (2)");
  strings.emplace_back("   {Endif}");
  strings.emplace_back("           {Ifdef(ifno)}");
  strings.emplace_back("This line should not be echoed. (3)");
  strings.emplace_back(" {Else}");
  strings.emplace_back("This line should be echoed (4)");
  strings.emplace_back("   {Endif}");
  strings.emplace_back("");
  strings.emplace_back("$ Test ifndef - else lines");
  strings.emplace_back(" {Ifndef(ifyes)}");
  strings.emplace_back("This line should not be echoed. (5)");
  strings.emplace_back("  {Else}");
  strings.emplace_back("This line should be echoed (6)");
  strings.emplace_back("   {Endif}");
  strings.emplace_back("    {Ifndef(ifno)}");
  strings.emplace_back("This line should be echoed. (7)");
  strings.emplace_back(" {Else}");
  strings.emplace_back("This line should not be echoed (8)");
  strings.emplace_back("  {Endif}");
  strings.emplace_back("$ Lines a, b, c, d, 1, 4, 6, 7 should be echoed");
  strings.emplace_back("$ Check line counting -- should be on line 78: {Parse Error}");
  strings.emplace_back("{ifdef(ifyes)} {This should be an error}");
  strings.emplace_back("{endif}");
  strings.emplace_back("");
  strings.emplace_back(
      "$ ========================================================================");
  strings.emplace_back("$ Test if lines");
  strings.emplace_back("{if(sqrt(4) == 2)}");
  strings.emplace_back("  This line should be echoed. (a)");
  strings.emplace_back("{endif}");
  strings.emplace_back("  This line should be echoed. (b)");
  strings.emplace_back("{if(sqrt(2) == 2 || sqrt(3) == 2)}");
  strings.emplace_back("This line should not be echoed");
  strings.emplace_back("{endif}");
  strings.emplace_back("This line should be echoed. (c)");
  strings.emplace_back("");
  strings.emplace_back("{diff = sqrt(3)*sqrt(3) - 3}");
  strings.emplace_back("$ Test if - else lines");
  strings.emplace_back("{if(sqrt(3)*sqrt(3) - 3 == diff)}");
  strings.emplace_back(" complex if works");
  strings.emplace_back("{else}");
  strings.emplace_back(" complex if does not work");
  strings.emplace_back("{endif}");
  strings.emplace_back("");
  strings.emplace_back("{if (sqrt(4) == 2)}");
  strings.emplace_back(" {if (sqrt(9) == 3)}");
  strings.emplace_back("  {if (sqrt(16) == 4)}");
  strings.emplace_back("    square roots work");
  strings.emplace_back("  {else}");
  strings.emplace_back("    sqrt(16) does not work");
  strings.emplace_back("  {endif}");
  strings.emplace_back(" {else}");
  strings.emplace_back("   sqrt(9) does not work");
  strings.emplace_back(" {endif}");
  strings.emplace_back("{else}");
  strings.emplace_back("  sqrt(4) does not work");
  strings.emplace_back("{endif}");
  strings.emplace_back("");
  strings.emplace_back("{v1 = 1} {v2 = 2}");
  strings.emplace_back("{if (v1 == v2)}");
  strings.emplace_back("  Bad if");
  strings.emplace_back("  {if (v1 != v2)}");
  strings.emplace_back("   should not see (1)");
  strings.emplace_back("  {else}");
  strings.emplace_back("   should not see (2)");
  strings.emplace_back("  {endif}");
  strings.emplace_back("   should not see (3)");
  strings.emplace_back("{else}");
  strings.emplace_back("  {if (v1 != v2)}");
  strings.emplace_back("   good nested if");
  strings.emplace_back("  {else}");
  strings.emplace_back("   bad nested if");
  strings.emplace_back("  {endif}");
  strings.emplace_back("  good");
  strings.emplace_back("  make sure it is still good");
  strings.emplace_back("{endif}");
  strings.emplace_back(
      "$ ========================================================================");
  strings.emplace_back("$ Test switch");
  strings.emplace_back("{switch(PI)}");
  strings.emplace_back("This is in a switch, but prior to any case, it should not run");
  strings.emplace_back("{a = 0.5} Should not execute");
  strings.emplace_back("");
  strings.emplace_back("{case (1)}");
  strings.emplace_back("This should not echo");
  strings.emplace_back("{a = 1}");
  strings.emplace_back("");
  strings.emplace_back("{case (2)}");
  strings.emplace_back("This should not echo");
  strings.emplace_back("{a = 2}");
  strings.emplace_back("");
  strings.emplace_back("{case (PI)}");
  strings.emplace_back("This should echo");
  strings.emplace_back("{a = PI}");
  strings.emplace_back("");
  strings.emplace_back("{case (PI)}");
  strings.emplace_back("This should not echo since a previous case matched.");
  strings.emplace_back("{a = 2}");
  strings.emplace_back("");
  strings.emplace_back("{default}");
  strings.emplace_back("{a=4}");
  strings.emplace_back("");
  strings.emplace_back("{endswitch}");
  strings.emplace_back("");
  strings.emplace_back("This should be equal to PI --  {PI}");
  strings.emplace_back("$ Test int and [] (shortcut for int)");
  strings.emplace_back("{int(5.01)}        {int(-5.01)}");
  strings.emplace_back("{[5.01]}   {[-5.01]}");
  strings.emplace_back("");
  strings.emplace_back("$ Test looping - print sin, cos from 0 to 90 by 5");
  strings.emplace_back("{_angle = -5}");
  strings.emplace_back("{Loop(19)}");
  strings.emplace_back(
      "{_angle += 5}    {_sa=sind(_angle)}      {_ca=cosd(_angle)} {hypot(_sa, _ca)} ");
  strings.emplace_back("{EndLoop}");
  strings.emplace_back("");
  strings.emplace_back("$$$$ Test formatting and string concatenation");
  strings.emplace_back("{_i = 0} {_SAVE = _FORMAT}");
  strings.emplace_back("{loop(20)}");
  strings.emplace_back(
      R("{IO(++_i)} Using the format {_FORMAT = " %." // tostring(_i) // " g "}, PI = {PI}"));
  strings.emplace_back("{endloop}");
  strings.emplace_back("Reset format to default: {_FORMAT = _SAVE}");
  strings.emplace_back("");
  strings.emplace_back("$$$$ Test string rescanning and executing");
  strings.emplace_back("{ECHO(OFF)}");
  strings.emplace_back("{Test = '  This is line 1: {a = atan2(0,-1)}");
  strings.emplace_back("        This is line 2: {sin(a/4)}");
  strings.emplace_back("   This is line 3: {cos(a/4)}'}");
  strings.emplace_back(R("{Test2 = 'This has an embedded string: {T = " This is a string "}'}"));
  strings.emplace_back("{ECHO(ON)}");
  strings.emplace_back("Original String:");
  strings.emplace_back("{Test}");
  strings.emplace_back("Rescanned String:");
  strings.emplace_back("{rescan(Test)} ");
  strings.emplace_back("Original String:");
  strings.emplace_back("{Test2}");
  strings.emplace_back("Print Value of variable T = {T}");
  strings.emplace_back("Rescanned String:");
  strings.emplace_back("{rescan(Test2)} ");
  strings.emplace_back("Print Value of variable T = {T}");
  strings.emplace_back("");
  strings.emplace_back(R("Original String: {t1 = " atan2(0, -1) "}"));
  strings.emplace_back("Executed String: {execute(t1)}");
  strings.emplace_back("");
  strings.emplace_back(R("string = {_string = " one two, three "}"));
  strings.emplace_back(R("delimiter " { _delm = " ," } ""));
  strings.emplace_back("word count = {word_count(_string,_delm)}");
  strings.emplace_back(R("second word = " { get_word(2, _string, _delm) } ""));
  strings.emplace_back("");
  strings.emplace_back(R("string = {_string = " (one two, three * four - five"}"));
  strings.emplace_back(R("delimiter " {
    _delm = " ,(*-" } ""));
  strings.emplace_back("word count = {word_count(_string,_delm)}");
       strings.emplace_back(R("second word = "{
    get_word(2, _string, _delm)}""));
  strings.emplace_back("");
  strings.emplace_back("");
  strings.emplace_back(R("string = {_string = " one two, three"}"));
  strings.emplace_back(R("delimiter " {
    _delm = " ," } ""));
  strings.emplace_back("word count = { iwords = word_count(_string,_delm)}");
  strings.emplace_back("");
  strings.emplace_back("{_n = 0}");
  strings.emplace_back("{loop(iwords)}");
  strings.emplace_back(R("word {++_n} = "{
    get_word(_n, _string, _delm)}""));
  strings.emplace_back("{endloop}");
  strings.emplace_back("");
  strings.emplace_back("$ Check parsing of escaped braces...");
  strings.emplace_back("\\{ int a = b + {PI/2} \\}");
  strings.emplace_back("\\{ \\}");
  strings.emplace_back(" ");

  bool result = aprepro.parse_strings(strings, "My list of strings");

  if (result) {
    std::string res_str = aprepro.parsing_results().str();
    std::cout << res_str;
  }

  aprepro.clear_results();
}
