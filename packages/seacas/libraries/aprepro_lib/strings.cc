// Copyright(C) 1999-2020 National Technology & Engineering Solutions
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
  strings.push_back("$ Test program for Aprepro");
  strings.push_back("$");
  strings.push_back("Test number representations");
  strings.push_back("{1}        {10e-1} {10.e-1}        {.1e+1} {.1e1}");
  strings.push_back("{1}        {10E-1} {10.E-1}        {.1E+1} {.1E1}");
  strings.push_back("   ");
  strings.push_back("Test assign statements:");
  strings.push_back("{_a = 5}   {b=_a}  $ Should print 5 5");
  strings.push_back("{_a +=b}   {_a}    $ Should print 10 10");
  strings.push_back("{_a -=b}   {_a}    $ Should print 5 5");
  strings.push_back("{_a *=b}   {_a}    $ Should print 25 25");
  strings.push_back("{_a /=b}   {_a}    $ Should print 5 5");
  strings.push_back("{_a ^=b}   {_a}    $ Should print 3125 3125");
  strings.push_back("{_a = b}   {_a**=b}        {_a}    $ Should print 5 3125 3125");
  strings.push_back("");
  strings.push_back("Test trigonometric functions (radians)");
  strings.push_back("{pi = d2r(180)} {atan2(0,-1)} {4*atan(1.0)} $ Three values of pi");
  strings.push_back("{_a = sin(pi/4)}   {pi-4*asin(_a)} $ sin(pi/4)");
  strings.push_back("{_b = cos(pi/4)}   {pi-4*acos(_b)} $ cos(pi/4)");
  strings.push_back("{_c = tan(pi/4)}   {pi-4*atan(_c)} $ tan(pi/4)");
  strings.push_back("");
  strings.push_back("Test trigonometric functions (degrees)");
  strings.push_back("{r2d(pi)}  {pid = atan2d(0,-1)}    {4 * atand(1.0)}");
  strings.push_back("{ad = sind(180/4)} {180-4*asind(ad)}       $ sin(180/4)");
  strings.push_back("{bd = cosd(180/4)} {180-4*acosd(bd)}       $ cos(180/4)");
  strings.push_back("{cd = tand(180/4)} {180-4*atand(cd)}       $ tan(180/4)");
  strings.push_back("");
  strings.push_back("Test max, min, sign, dim, abs");
  strings.push_back("{pmin = min(0.5, 1.0)}     {nmin = min(-0.5, -1.0)} $ Should be 0.5, -1");
  strings.push_back("{pmax = max(0.5, 1.0)}     {nmax = max(-0.5, -1.0)} $ Should be 1.0, -0.5");
  strings.push_back("{zero = 0} {sign(0.5, zero) + sign(0.5, -zero)}    $ Should be 0 1");
  strings.push_back("{nonzero = 1} {sign(0.5, nonzero) + sign(0.5, -nonzero)} $ Should be 1 0");
  strings.push_back("{dim(5.5, 4.5)}    {dim(4.5, 5.5)} $ Should be 1 0");
  strings.push_back("");
  strings.push_back("{ifyes = 1} {ifno = 0}");
  strings.push_back("$ Test ternary...");
  strings.push_back("{ifyes == 1 ? \"correct\" : \"incorrect\"}");
  strings.push_back("{ifno == 1 ? \"incorrect\" : \"correct\"}");
  strings.push_back("");
  strings.push_back("$ Test ifdef lines");
  strings.push_back("   {Ifdef(ifyes)}");
  strings.push_back("This line should be echoed. (a)");
  strings.push_back(" {Endif}");
  strings.push_back("This line should be echoed. (b)");
  strings.push_back("     {Ifdef(ifno)}");
  strings.push_back("This line should not be echoed");
  strings.push_back("    {Endif}");
  strings.push_back("This line should be echoed. (c)");
  strings.push_back("  {Ifdef(ifundefined)}");
  strings.push_back("This line should not be echoed");
  strings.push_back("        {Endif}");
  strings.push_back("This line should be echoed. (d)");
  strings.push_back("");
  strings.push_back("$ Test ifdef - else lines");
  strings.push_back("             {Ifdef(ifyes)}");
  strings.push_back("This line should be echoed. (1)");
  strings.push_back("                   {Else}");
  strings.push_back("This line should not be echoed (2)");
  strings.push_back("   {Endif}");
  strings.push_back("           {Ifdef(ifno)}");
  strings.push_back("This line should not be echoed. (3)");
  strings.push_back(" {Else}");
  strings.push_back("This line should be echoed (4)");
  strings.push_back("   {Endif}");
  strings.push_back("");
  strings.push_back("$ Test ifndef - else lines");
  strings.push_back(" {Ifndef(ifyes)}");
  strings.push_back("This line should not be echoed. (5)");
  strings.push_back("  {Else}");
  strings.push_back("This line should be echoed (6)");
  strings.push_back("   {Endif}");
  strings.push_back("    {Ifndef(ifno)}");
  strings.push_back("This line should be echoed. (7)");
  strings.push_back(" {Else}");
  strings.push_back("This line should not be echoed (8)");
  strings.push_back("  {Endif}");
  strings.push_back("$ Lines a, b, c, d, 1, 4, 6, 7 should be echoed");
  strings.push_back("$ Check line counting -- should be on line 78: {Parse Error}");
  strings.push_back("{ifdef(ifyes)} {This should be an error}");
  strings.push_back("{endif}");
  strings.push_back("");
  strings.push_back("$ ========================================================================");
  strings.push_back("$ Test if lines");
  strings.push_back("{if(sqrt(4) == 2)}");
  strings.push_back("  This line should be echoed. (a)");
  strings.push_back("{endif}");
  strings.push_back("  This line should be echoed. (b)");
  strings.push_back("{if(sqrt(2) == 2 || sqrt(3) == 2)}");
  strings.push_back("This line should not be echoed");
  strings.push_back("{endif}");
  strings.push_back("This line should be echoed. (c)");
  strings.push_back("");
  strings.push_back("{diff = sqrt(3)*sqrt(3) - 3}");
  strings.push_back("$ Test if - else lines");
  strings.push_back("{if(sqrt(3)*sqrt(3) - 3 == diff)}");
  strings.push_back(" complex if works");
  strings.push_back("{else}");
  strings.push_back(" complex if does not work");
  strings.push_back("{endif}");
  strings.push_back("");
  strings.push_back("{if (sqrt(4) == 2)}");
  strings.push_back(" {if (sqrt(9) == 3)}");
  strings.push_back("  {if (sqrt(16) == 4)}");
  strings.push_back("    square roots work");
  strings.push_back("  {else}");
  strings.push_back("    sqrt(16) does not work");
  strings.push_back("  {endif}");
  strings.push_back(" {else}");
  strings.push_back("   sqrt(9) does not work");
  strings.push_back(" {endif}");
  strings.push_back("{else}");
  strings.push_back("  sqrt(4) does not work");
  strings.push_back("{endif}");
  strings.push_back("");
  strings.push_back("{v1 = 1} {v2 = 2}");
  strings.push_back("{if (v1 == v2)}");
  strings.push_back("  Bad if");
  strings.push_back("  {if (v1 != v2)}");
  strings.push_back("   should not see (1)");
  strings.push_back("  {else}");
  strings.push_back("   should not see (2)");
  strings.push_back("  {endif}");
  strings.push_back("   should not see (3)");
  strings.push_back("{else}");
  strings.push_back("  {if (v1 != v2)}");
  strings.push_back("   good nested if");
  strings.push_back("  {else}");
  strings.push_back("   bad nested if");
  strings.push_back("  {endif}");
  strings.push_back("  good");
  strings.push_back("  make sure it is still good");
  strings.push_back("{endif}");
  strings.push_back("$ ========================================================================");
  strings.push_back("$ Test switch");
  strings.push_back("{switch(PI)}");
  strings.push_back("This is in a switch, but prior to any case, it should not run");
  strings.push_back("{a = 0.5} Should not execute");
  strings.push_back("");
  strings.push_back("{case (1)}");
  strings.push_back("This should not echo");
  strings.push_back("{a = 1}");
  strings.push_back("");
  strings.push_back("{case (2)}");
  strings.push_back("This should not echo");
  strings.push_back("{a = 2}");
  strings.push_back("");
  strings.push_back("{case (PI)}");
  strings.push_back("This should echo");
  strings.push_back("{a = PI}");
  strings.push_back("");
  strings.push_back("{case (PI)}");
  strings.push_back("This should not echo since a previous case matched.");
  strings.push_back("{a = 2}");
  strings.push_back("");
  strings.push_back("{default}");
  strings.push_back("{a=4}");
  strings.push_back("");
  strings.push_back("{endswitch}");
  strings.push_back("");
  strings.push_back("This should be equal to PI --  {PI}");
  strings.push_back("$ Test int and [] (shortcut for int)");
  strings.push_back("{int(5.01)}        {int(-5.01)}");
  strings.push_back("{[5.01]}   {[-5.01]}");
  strings.push_back("");
  strings.push_back("$ Test looping - print sin, cos from 0 to 90 by 5");
  strings.push_back("{_angle = -5}");
  strings.push_back("{Loop(19)}");
  strings.push_back(
      "{_angle += 5}    {_sa=sind(_angle)}      {_ca=cosd(_angle)} {hypot(_sa, _ca)} ");
  strings.push_back("{EndLoop}");
  strings.push_back("");
  strings.push_back("$$$$ Test formatting and string concatenation");
  strings.push_back("{_i = 0} {_SAVE = _FORMAT}");
  strings.push_back("{loop(20)}");
  strings.push_back(
      "{IO(++_i)} Using the format {_FORMAT = \"%.\" // tostring(_i) // \"g\"}, PI = {PI}");
  strings.push_back("{endloop}");
  strings.push_back("Reset format to default: {_FORMAT = _SAVE}");
  strings.push_back("");
  strings.push_back("$$$$ Test string rescanning and executing");
  strings.push_back("{ECHO(OFF)}");
  strings.push_back("{Test = '  This is line 1: {a = atan2(0,-1)}");
  strings.push_back("        This is line 2: {sin(a/4)}");
  strings.push_back("   This is line 3: {cos(a/4)}'}");
  strings.push_back("{Test2 = 'This has an embedded string: {T = \"This is a string\"}'}");
  strings.push_back("{ECHO(ON)}");
  strings.push_back("Original String:");
  strings.push_back("{Test}");
  strings.push_back("Rescanned String:");
  strings.push_back("{rescan(Test)} ");
  strings.push_back("Original String:");
  strings.push_back("{Test2}");
  strings.push_back("Print Value of variable T = {T}");
  strings.push_back("Rescanned String:");
  strings.push_back("{rescan(Test2)} ");
  strings.push_back("Print Value of variable T = {T}");
  strings.push_back("");
  strings.push_back("Original String: {t1 = \"atan2(0,-1)\"}");
  strings.push_back("Executed String: {execute(t1)}");
  strings.push_back("");
  strings.push_back("string = {_string = \" one two, three\"}");
  strings.push_back("delimiter \"{_delm = \" ,\"}\"");
  strings.push_back("word count = {word_count(_string,_delm)}");
  strings.push_back("second word = \"{get_word(2,_string,_delm)}\"");
  strings.push_back("");
  strings.push_back("string = {_string = \" (one two, three * four - five\"}");
  strings.push_back("delimiter \"{_delm = \" ,(*-\"}\"");
  strings.push_back("word count = {word_count(_string,_delm)}");
  strings.push_back("second word = \"{get_word(2,_string,_delm)}\"");
  strings.push_back("");
  strings.push_back("");
  strings.push_back("string = {_string = \" one two, three\"}");
  strings.push_back("delimiter \"{_delm = \" ,\"}\"");
  strings.push_back("word count = { iwords = word_count(_string,_delm)}");
  strings.push_back("");
  strings.push_back("{_n = 0}");
  strings.push_back("{loop(iwords)}");
  strings.push_back("word {++_n} = \"{get_word(_n,_string,_delm)}\"");
  strings.push_back("{endloop}");
  strings.push_back("");
  strings.push_back("$ Check parsing of escaped braces...");
  strings.push_back("\\{ int a = b + {PI/2} \\}");
  strings.push_back("\\{ \\}");
  strings.push_back(" ");

  bool result = aprepro.parse_strings(strings, "My list of strings");

  if (result) {
    std::string res_str = aprepro.parsing_results().str();
    std::cout << res_str;
  }

  aprepro.clear_results();
}
