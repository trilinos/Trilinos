#!/usr/bin/env python

"""
Generate a .csv file for import into Drupal database supporting the Trilinos web
site.  Data is read from a cmake-generated XML file and optionally from a Drupal-
feeds-generated .csv file.
"""

# System imports
from optparse import *
from xml.dom.minidom import parse
import xml.dom.minidom
import csv
import re
import sys
import os

###############################################################################
def extract_version(filepath):
    """
    Read the Trilinos version file, parse it, and extract the version.  In
    special cases, override this version string with "Unspec" (for unspecified),
    and "Devel" for the development branch.  Return the version string.
    """
    # Read the version file, parse it and extract the version
    version = "Unspec"
    for line in open(filepath,'r').readlines():
        if line.startswith("SET"):
            if line.find("Trilinos_VERSION_STRING") > 0:
                version = line.split('"')[1]
    if version.find("Dev") > 0: version = "Devel"
    return version

###############################################################################
def package_id(package, version):
    """
    Return a unique package ID, given the package name and version string.
    Currently, this is a concatenation of name, "-", and version.
    """
    return package + "-" + version

###############################################################################
def parse_xml_file(filepath, version):
    """
    Parse the Trilinos package dependencies XML file and return a database of
    software packages: both Trilinos packages and third-party libraries used by
    Trilinos.  The DB is a python dictionary with package IDs for keys and
    "package dictionaries" for values.  A package ID is a unique string based on
    the package name and version number, obtained from the function
    package_id(name,version).  A package dictionary is a python dictionary that
    is guaranteed to contain the following keys: name, version, parent,
    sub_packages, required_dep, optional_dep.  The name, version and parent
    values are strings; the remaining values are lists of strings.

    This function conducts only a first pass through the data.  The parent
    fields will be accurate package IDs, but the sub_package fields will be
    empty.  Also, the required_dep and optional_dep lists will contain package
    names, not package IDs.  This data must be corrected on subsequent passes
    through the data with different functions.
    """

    # Set up the Trilinos package dictionary and insert it as the first entry
    # for the package database
    trilinos = u"Trilinos"
    trilinosID = package_id(trilinos, version)
    package_dict = {"name"         : trilinos,
                    "version"      : version,
                    "parent"       : "",
                    "sub_packages" : [],
                    "required_dep" : [],
                    "optional_dep" : []
                    }
    package_db = {trilinosID : package_dict}

    # Define th expected attributes in the XML file
    xml_attributes = ["LIB_REQUIRED_DEP_PACKAGES",
                      "TEST_REQUIRED_DEP_PACKAGES",
                      "LIB_REQUIRED_DEP_TPLS",
                      "TEST_REQUIRED_DEP_TPLS",
                      "LIB_OPTIONAL_DEP_PACKAGES",
                      "TEST_OPTIONAL_DEP_PACKAGES",
                      "LIB_OPTIONAL_DEP_TPLS",
                      "TEST_OPTIONAL_DEP_TPLS",
                      "ParentPackage"
                      ]
 
    DOMTree = xml.dom.minidom.parse(filepath)
    collection = DOMTree.documentElement

    if collection.hasAttribute("project"):
        project = collection.getAttribute("project")
        if project != trilinos:
            raise ValueError("XML file does not specify Trilinos dependencies")

    # Get all the packages in the collection
    packages = collection.getElementsByTagName("Package")
    
    # For each package, generate a dictionary that contains all of the build
    # information provided by the dependencies file, and then insert that
    # dictionary into the package database
    for package in packages:

        # Initialize the package dictionary, using the Trilinos version as the
        # default version and the Trilinos ID as the default parent
        package_dict = {"name"         : "",
                        "version"      : version,
                        "parent"       : trilinosID,
                        "sub_packages" : [],
                        "required_dep" : [],
                        "optional_dep" : []
                        }

        # Package name
        if package.hasAttribute("name"):
            package_name = package.getAttribute("name")
            package_dict["name"] = package_name
        else:
            print "No package name found"
            next

        # Initialize the dependencies and sub-packages
        package_dict["sub_packages"] = []
        package_dict["required_dep"] = [] 
        package_dict["optional_dep"] = []

        # Loop over the valid attributes
        for attribute in xml_attributes:
            element = package.getElementsByTagName(attribute)[0]
            if element:
                value = element.getAttribute("value")
                if value:
                    packages = value.split(",")
                    if "REQUIRED" in attribute:
                        package_dict["required_dep"].extend(packages)
                    elif "OPTIONAL" in attribute:
                        package_dict["optional_dep"].extend(packages)
                    elif attribute == "ParentPackage":
                        package_dict["parent"] = package_id(value, version)
        
        package_db[package_id(package_name, version)] = package_dict

    return package_db

###############################################################################
def get_trilinos_id(package_db):
    """
    Given the software package database, find the key that corresponds to the
    Trilinos ID.  This will be a unique string that contains both "Trilinos" and
    the version string.  Return this package ID.
    """
    trilinos = u"Trilinos"
    trilinosID = None
    for package in package_db.keys():
        if package.startswith(trilinos):
            trilinosID = package
            break
    if trilinosID is None:
        raise ValueError("Trilinos not found in database")
    return trilinosID

###############################################################################
def get_trilinos_version(package_db):
    """
    Given the software package database, find the key that corresponds to the
    Trilinos ID.  Extract from this the Trilinos version string.  Return this
    version string.
    """
    trilinosID = get_trilinos_id(package_db)
    return trilinosID.split("-")[1]

###############################################################################
def get_subpackages(package_db):
    """
    Given the package database produced by parse_xml_file(), translate the
    parent information into accurate sub-package entries.

    This function should be the second pass through the package data when
    forming the package database.
    """
    for packageID, packageData in package_db.items():
        # The parent ID is guaranteed to exist
        parentID = packageData["parent"]
        # If the parent ID is the null string, then the package being considered
        # is the top-level Trilinos package.  Otherwise, the parent ID should be
        # either the Trilinos ID, or an ID of a valid Trilinos package.
        if parentID:
            package_db[parentID]["sub_packages"].append(packageID)

###############################################################################
def extract_tpls(package_db):
    """
    Given the package database, extract the list of third-party libraries.

    This is achieved by building two sets: the set of trilinos package names,
    and the set of all required and optional dependencies.  The TPLs are the set
    that remains when you subtract the trilinos packages from the dependencies.
    """
    # Get the set of Trilinos package names (not IDs) and the set of all
    # dependencies
    trilinos_package_names = set()
    dependencies           = set()
    for package in package_db.values():
        trilinos_package_names.add(package["name"])
        dependencies.update(package["required_dep"])
        dependencies.update(package["optional_dep"])

    # Return the sorted TPLs
    tpls = list(dependencies.difference(trilinos_package_names))
    tpls.sort()
    return tpls

###############################################################################
def process_dependencies(package_db):
    """
    This function should be the third pass through the package data when forming
    the package database.

    On entry, the package database should have only entries for Trilinos
    packages, and those packages' required and optional dependencies should be
    unordered lists of package names (not package IDs), and can contain
    duplicates.

    On exit, the package database will be updated so that Trilinos packages'
    required and optional dependencies are sorted lists of unique package IDs
    (not just names).  Also, the package database will be extended to include
    third-party libraries.

    Note that for now, TPL version strings will always be "Unspec".
    """
    trilinos_version = get_trilinos_version(package_db)
    tpl_version      = "Unspec"
    tpl_names        = extract_tpls(package_db)

    # Change dependency lists from package names to package IDs
    for packageData in package_db.values():
        for dep in ("required_dep", "optional_dep"):
            dependencies = packageData[dep]
            for i in range(len(dependencies)):
                if dependencies[i] in tpl_names:
                    dependencies[i] = package_id(dependencies[i],
                                                 tpl_version)
                else:
                    dependencies[i] = package_id(dependencies[i],
                                                 trilinos_version)
            # Remove duplicates and sort the dependencies
            new_dep = list(set(dependencies))
            new_dep.sort()
            packageData[dep] = new_dep

    # Add TPLs to the package database
    for tpl_name in tpl_names:
        tpl_dict = {"name"         : tpl_name,
                    "version"      : tpl_version,
                    "parent"       : "",
                    "sub_packages" : [],
                    "required_dep" : [],
                    "optional_dep" : []
                    }
        tpl_id = package_id(tpl_name, tpl_version)
        package_db[tpl_id] = tpl_dict

###############################################################################
def process_xml_file(filepath, version):
    """
    Given a path to the Trilinos package dependencies XML file and a Trilinos
    version string, return a complete package database.
    """
    package_db = parse_xml_file(filepath, version)
    get_subpackages(package_db)
    process_dependencies(package_db)
    return package_db

###############################################################################
def default_csv_output_file(version):
    return "Trilinos-" + version + ".csv"

###############################################################################
def stringify_package_dict(package_dict):
    """
    Return a new package dictionary equivalent to the given package dictionary
    in which all Unicode strings have been converted to simple strings and each
    list of (possibly Unicode) strings has been converted to a single, simple,
    comma-separated string.
    """
    new_dict = {"name"         : str(package_dict["name"   ]),
                "version"      : str(package_dict["version"]),
                "parent"       : str(package_dict["parent" ]),
                "sub_packages" : str(",".join(package_dict["sub_packages"])),
                "required_dep" : str(",".join(package_dict["required_dep"])),
                "optional_dep" : str(",".join(package_dict["optional_dep"]))
                }
    return new_dict

###############################################################################
def database_stats(package_db):
    """
    Return a dictionary storing statistical information about the given package
    database.
    """
    num_packages              = len(package_db)
    num_trilinos_packages     = 0
    num_trilinos_sub_packages = 0
    num_tpls                  = 0
    for packageData in package_db.values():
        name   = packageData["name"  ]
        parent = packageData["parent"]
        if name == "Trilinos":
            pass
        elif parent.startswith("Trilinos"):
            num_trilinos_packages += 1
        elif parent == "":
            num_tpls += 1
        else:
            num_trilinos_sub_packages += 1
    return {"Total number of packages"        : num_packages,
            "Number of Trilinos packages"     : num_trilinos_packages,
            "Number of Trilinos sub-packages" : num_trilinos_sub_packages,
            "Number of third-party libraries" : num_tpls
            }

###############################################################################
def write_to_csv(filepath, package_db):
    """
    Use the csv module to write the contents of the package database to a
    comma-separated value file.
    """
    # Setup
    trilinosID    = get_trilinos_id(package_db)
    trilinos_dict = package_db[trilinosID]
    fieldnames    = ["package_id",
                     "name",
                     "version",
                     "parent",
                     "sub_packages",
                     "required_dep",
                     "optional_dep"
                     ]

    # Write the file
    with open(filepath, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for packageID, package in package_db.items():
            row = stringify_package_dict(package)
            row["package_id"] = packageID
            writer.writerow(row)

###############################################################################
def main(trilinos_dependencies,
         trilinos_version_file,
         csv_input_file,
         csv_output_file,
         verbose):

    # Print info
    if verbose:
        print "Trilinos dependencies file:", trilinos_dependencies
        print "Trilinos version file     :", trilinos_version_file

    # Extract the Trilinos version number
    version = extract_version(trilinos_version_file)
    if verbose: print "Trilinos version          : '%s'" % version

    # Process the XML file to generate a package database
    package_db = process_xml_file(trilinos_dependencies, version)

    # Report the statistics
    if verbose > 1:
        stats = database_stats(package_db)
        headings = stats.keys()
        headings.sort()
        lengths = [len(x) for x in headings]
        lengths.sort()
        format = r"    %-" + str(lengths[-1]) + r"s: %5d"
        for heading in headings:
            print format % (heading, stats[heading])

    # Write the CSV output file
    if csv_output_file is None:
        csv_output_file = default_csv_output_file(version)
    if verbose: print "CSV output file: '%s'" % csv_output_file,
    write_to_csv(csv_output_file, package_db)
    if verbose: print "... written"

###############################################################################
if (__name__ == "__main__"):

    # Parse the command-line arguments
    parser = OptionParser()
    parser.add_option("-i", "--csvinput", action="store", dest="csv_input",
                      default=None, help="CSV input file name")
    parser.add_option("-o", "--csvoutput", action="store", dest="csv_output",
                      default=None, help="CSV output file name")
    parser.add_option("-v", "--verbosity", type="int", dest="verbosity",
                      default=1, help="set the verbosity level [default 1]")
    options,args = parser.parse_args()

    # Get the Trilinos source directory
    if len(args) >= 1:
        trilinos_src = args[0]
    else:
        print "usage:", sys.argv[0], "trilinos_src"
        sys.exit(-1)

    # Construct the dependencies file name
    trilinos_dependencies = os.path.join(trilinos_src,
                                         "cmake",
                                         "dependencies",
                                         "TrilinosPackageDependencies.xml")

    # Construct the version file name
    trilinos_version_file = os.path.join(trilinos_src,
                                         "Version.cmake")

    # Get the CSV input and output files
    csv_input_file  = options.csv_input
    csv_output_file = options.csv_output

    # Check for file existence
    files_to_check = [trilinos_dependencies, trilinos_version_file]
    if csv_input_file:
        files_to_check.append(csv_input_file)
    for file in files_to_check:
        if not os.path.isfile(file):
            print "'%s' not found" % file
            sys.exit(-2)

    # Call the main routine
    main(trilinos_dependencies,
         trilinos_version_file,
         csv_input_file,
         csv_output_file,
         options.verbosity)
