#!/bin/env python


from TrilinosDependencies import getTrilinosDependenciesFromXmlFile

trilinosDependencies = getTrilinosDependenciesFromXmlFile(
  'data/TrilinosPackageDependencies.xml')

print "\ntrilinosDependencies:\n", trilinosDependencies
