Instructions for using CppUnitLite with MS VC++:

Extract all these filed to a directory (CppUnitLite)
Move StackMain.cpp, StackTest.cpp and Stack.h to some other directory
Create a workspace
Create a project for a static library called CppUnitLite, turn off any pre-compiled headers stuff
Add all the CppUnitLite cpp and h files to the CppUnitLite project
Get it to compile
Create another project (Stack)
Add StackMain.cpp, StackTest.cpp and Stack.h to the project
Make the Stack project depend on CppUnitLite
Build it
Run it
Evolve it

