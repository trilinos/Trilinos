#!/bin/env python


from TrilinosDependencies import getTrilinosDependenciesFromXmlFile

trilinosDependencies = getTrilinosDependenciesFromXmlFile()

#print "\ntrilinosDependencies:\n", trilinosDependencies

trilinosLibDependenciesTable = trilinosDependencies.createRawTable(True)

print "\ntrilinosLibDependenciesTable:\n", trilinosLibDependenciesTable

trilinosDependenciesTable = trilinosDependencies.createRawTable(False)

print "\ntrilinosDependenciesTable:\n", trilinosDependenciesTable
