// Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "aprepro.h" // for Aprepro
#include <iostream>  // for cout, ostream, etc
#include <string>    // for string, operator<<
#include <vector>    // for vector

std::vector<std::string> build_strings();

int main(int, char **)
{
  SEAMS::Aprepro aprepro;

  aprepro.ap_options.warning_msg = false;

  std::vector<std::string> strings = build_strings();
  bool                     result  = aprepro.parse_strings(strings, "My list of strings");
  if (result) {
    std::string res_str = aprepro.parsing_results().str();
    std::cout << res_str;
  }
  aprepro.clear_results();
}

std::vector<std::string> build_strings()
{
  std::vector<std::string> strings;

  strings.emplace_back(R"($ Test program for Aprepro)");
  strings.emplace_back(R"($)");
  strings.emplace_back(R"(Test number representations)");
  strings.emplace_back(R"({1}	{10e-1}	{10.e-1}	{.1e+1}	{.1e1})");
  strings.emplace_back(R"({1}	{10E-1}	{10.E-1}	{.1E+1}	{.1E1})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(Test assign statements:)");
  strings.emplace_back(R"({_a = 5}	{b=_a}	$ Should print 5 5)");
  strings.emplace_back(R"({_a +=b}	{_a} 	$ Should print 10 10)");
  strings.emplace_back(R"({_a -=b}	{_a}	$ Should print 5 5)");
  strings.emplace_back(R"({_a *=b}	{_a}	$ Should print 25 25)");
  strings.emplace_back(R"({_a /=b}	{_a}	$ Should print 5 5)");
  strings.emplace_back(R"({_a ^=b}	{_a}	$ Should print 3125 3125)");
  strings.emplace_back(R"({_a = b}	{_a**=b}	{_a}	$ Should print 5 3125 3125)");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(Test trigonometric functions (radians))");
  strings.emplace_back(R"({pi = d2r(180)} {atan2(0,-1)} {4*atan(1.0)} $ Three values of pi)");
  strings.emplace_back(R"({_a = sin(pi/4)}	{pi-4*asin(_a)}	$ sin(pi/4))");
  strings.emplace_back(R"({_b = cos(pi/4)}	{pi-4*acos(_b)}	$ cos(pi/4))");
  strings.emplace_back(R"({_c = tan(pi/4)}	{pi-4*atan(_c)}	$ tan(pi/4))");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(Test trigonometric functions (degrees))");
  strings.emplace_back(R"({r2d(pi)}	{pid = atan2d(0,-1)}	{4 * atand(1.0)})");
  strings.emplace_back(R"({ad = sind(180/4)}	{180-4*asind(ad)}	$ sin(180/4))");
  strings.emplace_back(R"({bd = cosd(180/4)}	{180-4*acosd(bd)}	$ cos(180/4))");
  strings.emplace_back(R"({cd = tand(180/4)}	{180-4*atand(cd)}	$ tan(180/4))");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(Test max, min, sign, dim, abs)");
  strings.emplace_back(R"({pmin = min(0.5, 1.0)}	{nmin = min(-0.5, -1.0)} $ Should be 0.5, -1)");
  strings.emplace_back(R"({pmax = max(0.5, 1.0)}	{nmax = max(-0.5, -1.0)} $ Should be 1.0, -0.5)");
  strings.emplace_back(R"({zero = 0} {sign(0.5, zero) + sign(0.5, -zero)}	$ Should be 0 1)");
  strings.emplace_back(
      R"({nonzero = 1} {sign(0.5, nonzero) + sign(0.5, -nonzero)} $ Should be 1 0)");
  strings.emplace_back(R"({dim(5.5, 4.5)}	{dim(4.5, 5.5)}	$ Should be 1 0)");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({ifyes = 1} {ifno = 0})");
  strings.emplace_back(R"($ Test ternary...)");
  strings.emplace_back(R"({ifyes == 1 ? "correct" : "incorrect"})");
  strings.emplace_back(R"({ifno == 1 ? "incorrect" : "correct"})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test ifdef lines)");
  strings.emplace_back(R"(	{Ifdef(ifyes)})");
  strings.emplace_back(R"(This line should be echoed. (a))");
  strings.emplace_back(R"( {Endif})");
  strings.emplace_back(R"(This line should be echoed. (b))");
  strings.emplace_back(R"(     {Ifdef(ifno)})");
  strings.emplace_back(R"(This line should not be echoed)");
  strings.emplace_back(R"( 	 {Endif})");
  strings.emplace_back(R"(This line should be echoed. (c))");
  strings.emplace_back(R"(  {Ifdef(ifundefined)})");
  strings.emplace_back(R"(This line should not be echoed)");
  strings.emplace_back(R"(        {Endif})");
  strings.emplace_back(R"(This line should be echoed. (d))");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test ifdef - else lines)");
  strings.emplace_back(R"(             {Ifdef(ifyes)})");
  strings.emplace_back(R"(This line should be echoed. (1))");
  strings.emplace_back(R"(			{Else})");
  strings.emplace_back(R"(This line should not be echoed (2))");
  strings.emplace_back(R"(	{Endif})");
  strings.emplace_back(R"(		{Ifdef(ifno)})");
  strings.emplace_back(R"(This line should not be echoed. (3))");
  strings.emplace_back(R"( {Else})");
  strings.emplace_back(R"(This line should be echoed (4))");
  strings.emplace_back(R"(   {Endif})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test ifndef - else lines)");
  strings.emplace_back(R"( {Ifndef(ifyes)})");
  strings.emplace_back(R"(This line should not be echoed. (5))");
  strings.emplace_back(R"(  {Else})");
  strings.emplace_back(R"(This line should be echoed (6))");
  strings.emplace_back(R"(   {Endif})");
  strings.emplace_back(R"(    {Ifndef(ifno)})");
  strings.emplace_back(R"(This line should be echoed. (7))");
  strings.emplace_back(R"( {Else})");
  strings.emplace_back(R"(This line should not be echoed (8))");
  strings.emplace_back(R"(  {Endif})");
  strings.emplace_back(R"($ Lines a, b, c, d, 1, 4, 6, 7 should be echoed)");
  strings.emplace_back(R"($ Check line counting -- should be on line 78: )");
  strings.emplace_back(R"( )");
  strings.emplace_back(
      R"($ ========================================================================)");
  strings.emplace_back(R"($ Test string if lines)");
  strings.emplace_back(R"({if("Greg")})");
  strings.emplace_back(R"( This line should be echoed ("greg"))");
  strings.emplace_back(R"({else})");
  strings.emplace_back(R"( This line should not be echoed ("greg"))");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({empty=""})");
  strings.emplace_back(R"({if(empty)})");
  strings.emplace_back(R"( This line should not be echoed (empty))");
  strings.emplace_back(R"({else})");
  strings.emplace_back(R"( This line should be echoed (empty))");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"()");
  strings.emplace_back(
      R"($ ========================================================================)");
  strings.emplace_back(R"($ Test if lines)");
  strings.emplace_back(R"({if(sqrt(4) == 2)})");
  strings.emplace_back(R"(  This line should be echoed. (a))");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"(  This line should be echoed. (b))");
  strings.emplace_back(R"({if(sqrt(2) == 2 || sqrt(3) == 2)})");
  strings.emplace_back(R"( This line should not be echoed)");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"(This line should be echoed. (c))");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({diff = sqrt(3)*sqrt(3) - 3})");
  strings.emplace_back(R"($ Test if - else lines)");
  strings.emplace_back(R"({if(sqrt(3)*sqrt(3) - 3 == diff)})");
  strings.emplace_back(R"( complex if works)");
  strings.emplace_back(R"({else})");
  strings.emplace_back(R"( complex if does not work)");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({if (sqrt(4) == 2)})");
  strings.emplace_back(R"( {if (sqrt(9) == 3)})");
  strings.emplace_back(R"(  {if (sqrt(16) == 4)})");
  strings.emplace_back(R"(    square roots work)");
  strings.emplace_back(R"(  {else})");
  strings.emplace_back(R"(    sqrt(16) does not work)");
  strings.emplace_back(R"(  {endif})");
  strings.emplace_back(R"( {else})");
  strings.emplace_back(R"(   sqrt(9) does not work)");
  strings.emplace_back(R"( {endif})");
  strings.emplace_back(R"({else})");
  strings.emplace_back(R"(  sqrt(4) does not work)");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({v1 = 1} {v2 = 2})");
  strings.emplace_back(R"({if (v1 == v2)})");
  strings.emplace_back(R"(  Bad if)");
  strings.emplace_back(R"(  {if (v1 != v2)})");
  strings.emplace_back(R"(   should not see (1))");
  strings.emplace_back(R"(  {else})");
  strings.emplace_back(R"(   should not see (2))");
  strings.emplace_back(R"(  {endif})");
  strings.emplace_back(R"(   should not see (3))");
  strings.emplace_back(R"({else})");
  strings.emplace_back(R"(  {if (v1 != v2)})");
  strings.emplace_back(R"(   good nested if)");
  strings.emplace_back(R"(  {else})");
  strings.emplace_back(R"(   bad nested if)");
  strings.emplace_back(R"(  {endif})");
  strings.emplace_back(R"(  good)");
  strings.emplace_back(R"(  make sure it is still good)");
  strings.emplace_back(R"({endif})");
  strings.emplace_back(
      R"($ ========================================================================)");
  strings.emplace_back(R"($ Test switch)");
  strings.emplace_back(R"({switch(PI)})");
  strings.emplace_back(R"(This is in a switch, but prior to any case, it should not run)");
  strings.emplace_back(R"({a = 0.5} Should not execute)");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({case (1)})");
  strings.emplace_back(R"(This should not echo)");
  strings.emplace_back(R"({a = 1})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({case (2)})");
  strings.emplace_back(R"(This should not echo)");
  strings.emplace_back(R"({a = 2})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({case (PI)})");
  strings.emplace_back(R"(This should echo)");
  strings.emplace_back(R"({a = PI})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({case (PI)})");
  strings.emplace_back(R"(This should not echo since a previous case matched.)");
  strings.emplace_back(R"({a = 2})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({default})");
  strings.emplace_back(R"({a=4})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"({endswitch})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(This should be equal to PI --  {PI})");
  strings.emplace_back(R"($ Test int and [] (shortcut for int))");
  strings.emplace_back(R"({int(5.01)}	{int(-5.01)})");
  strings.emplace_back(R"({[5.01]}	{[-5.01]})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test looping - print sin, cos from 0 to 90 by 5)");
  strings.emplace_back(R"({Loop(19, _angle, 0, 5)})");
  strings.emplace_back(
      R"({_angle}	{_sa=sind(_angle)}	{_ca=cosd(_angle)} {hypot(_sa, _ca)} )");
  strings.emplace_back(R"({EndLoop})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($$$$ Test formatting and string concatenation)");
  strings.emplace_back(R"({_SAVE = _FORMAT})");
  strings.emplace_back(R"({loop(20)})");
  strings.emplace_back(
      R"({IO(__loop_1+1)} Using the format {_FORMAT = "%." // tostring(__loop_1+1) // "g"},	PI = {PI})");
  strings.emplace_back(R"({endloop})");
  strings.emplace_back(R"(Reset format to default: {_FORMAT = _SAVE})");
  strings.emplace_back(R"()");
  strings.emplace_back(
      R"($$$$ Test formatting using the `format` function.  _FORMAT is not modified)");
  strings.emplace_back(R"({loop(20)}")");
  strings.emplace_back(
      R"({__loop_1+1} Using the format {_f = "%." // tostring(__loop_1+1) // "f"},	PI = {format(PI,_f)})");
  strings.emplace_back(R"({endloop})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($$$$ Test string rescanning and executing)");
  strings.emplace_back(R"({ECHO(OFF)})");
  strings.emplace_back(R"({Test = '	This is line 1: {a = atan2(0,-1)})");
  strings.emplace_back(R"(        This is line 2: {sin(a/4)})");
  strings.emplace_back(R"(	This is line 3: {cos(a/4)}'})");
  strings.emplace_back(R"({Test2 = 'This has an embedded string: {T = "This is a string"}'})");
  strings.emplace_back(R"({ECHO(ON)})");
  strings.emplace_back(R"(Original String:)");
  strings.emplace_back(R"({Test})");
  strings.emplace_back(R"(Rescanned String:)");
  strings.emplace_back(R"({rescan(Test)} )");
  strings.emplace_back(R"(Original String:)");
  strings.emplace_back(R"({Test2})");
  strings.emplace_back(R"(Print Value of variable T = {T})");
  strings.emplace_back(R"(Rescanned String:)");
  strings.emplace_back(R"({rescan(Test2)} )");
  strings.emplace_back(R"(Print Value of variable T = {T})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"d(Original String: {t1 = "atan2(0,-1)"})d");
  strings.emplace_back(R"(Executed String: {execute(t1)})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(string = {_string = " one two, three"})");
  strings.emplace_back(R"(delimiter "{_delm = " ,"}")");
  strings.emplace_back(R"(word count = {word_count(_string,_delm)})");
  strings.emplace_back(R"(second word = "{get_word(2,_string,_delm)}")");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(string = {_string = " (one two, three * four - five"})");
  strings.emplace_back(R"(delimiter "{_delm = " ,(*-"}")");
  strings.emplace_back(R"(word count = {word_count(_string,_delm)})");
  strings.emplace_back(R"(second word = "{get_word(2,_string,_delm)}")");
  strings.emplace_back(R"()");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(string = {_string = " one two, three"})");
  strings.emplace_back(R"(delimiter "{_delm = " ,"}")");
  strings.emplace_back(R"(word count = { iwords = word_count(_string,_delm)})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"(	{loop(iwords, _n, 1)})");
  strings.emplace_back(R"(word {_n} = "{get_word(_n,_string,_delm)}")");
  strings.emplace_back(R"(   {endloop})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Check parsing of escaped braces...)");
  strings.emplace_back(R"(\{ int a = b + {PI/2} \})");
  strings.emplace_back(R"(\{ \})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test variable deletion)");
  strings.emplace_back(R"({new_var = sqrt(2) * sqrt(3)})");
  strings.emplace_back(R"({new_var})");
  strings.emplace_back(R"({delete("new_var")})");
  strings.emplace_back(R"({new_var}  This should print warning about undefined variable)");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test extract)");
  strings.emplace_back(R"({ex_found = extract("The test string is found", "test", "")})");
  strings.emplace_back(R"({ex_null  = extract("The test string is not found", "xxxx", "yyyy")})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($ Test string tokenization optimization)");
  strings.emplace_back(R"({list1 ='51,52,53,54,61,62,63,64'})");
  strings.emplace_back(R"({list2 ='71,72,73,74,81,82,83,84'})");
  strings.emplace_back(R"({loop(8, _i, 1)})");
  strings.emplace_back(
      R"(Word {_i} of list1 and list2 are {get_word(_i,list1,',')} and {get_word(_i,list2,',')})");
  strings.emplace_back(R"({endloop})");
  strings.emplace_back(R"()");
  strings.emplace_back(R"($$$$ Test double brace echo off/on)");
  strings.emplace_back(R"(Nothing further on line: {{"not echoed"}})");
  strings.emplace_back(
      R"(Noecho followed by non-parsing output: {{"not echoed"}}This should be echoed)");
  strings.emplace_back(
      R"(Echo, noecho setting variable, then echo that variable: {e="echo"}+{{d="echo"}}+{d})");
  strings.emplace_back(R"($End of test file)");

  return strings;
}
