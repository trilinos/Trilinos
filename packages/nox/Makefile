########################################################################
### Makefile
########################################################################

########################################################################
### Definitions 
########################################################################

CXXFLAGS = -g -I$(TRILINOS_HOME)/src/petra

########################################################################
### Executable Rules 
########################################################################

OBJS = NLS_Newton.o NLS_PetraVector.o NLS_Method.o NLS_ParameterList.o 

tester: NLS_ParameterList.o NLS_Parameter.o tester.o
	 $(CXX) $(CXXFLAGS) tester.o NLS_ParameterList.o NLS_Parameter.o -o tester

########################################################################
### Object File Rules
########################################################################

NLS_Newton.o : NLS_Newton.C NLS_Method.H NLS_Newton.H 
	$(CXX) -c $(CXXFLAGS) NLS_Newton.C -o $@

NLS_ParameterList.o : NLS_ParameterList.C NLS_Parameter.H NLS_ParameterList.H 
	$(CXX) -c $(CXXFLAGS) NLS_ParameterList.C -o $@

NLS_PetraVector.o : NLS_PetraVector.C NLS_PetraVector.H 
	$(CXX) -c $(CXXFLAGS) NLS_PetraVector.C -o $@

tester.o : tester.C NLS_Parameter.H NLS_ParameterList.H 
	$(CXX) -c $(CXXFLAGS) tester.C -o $@

NLS_GenericMethod.o : NLS_GenericMethod.C NLS_Function.H NLS_Group.H \
	NLS_Method.H NLS_Newton.H NLS_ParameterList.H 
	$(CXX) -c $(CXXFLAGS) NLS_GenericMethod.C -o $@

NLS_Parameter.o : NLS_Parameter.C 
	$(CXX) -c $(CXXFLAGS) NLS_Parameter.C -o $@

########################################################################
### Utilities
########################################################################



