#!/usr/bin/env python
import xml.etree.ElementTree as ET
import argparse


parser = argparse.ArgumentParser(description="Convert an extended xml input file to use Kokkos")
parser.add_argument("-o", "--output", dest="output",
                    default="out.xml", help="xml output file")
parser.add_argument('xml', metavar='xml', type=str, nargs='+',
                    help='an integer for the accumulator')

options = parser.parse_args()

if len(options.xml) > 1:
    outfiles = options.xml
else:
    outfiles = [options.output]

for fn, out in zip(options.xml, outfiles):
    tree = ET.parse(fn)
    root = tree.getroot()

    factoryList = root.findall("*[@name='Factories']")
    if len(factoryList) == 1:
        # extended format
        factoryList = factoryList[0]
        factoriesToDelete = []
        for factoryParams in factoryList.findall("ParameterList"):
            factory = factoryParams.findall("*[@name='factory']")
            assert len(factory) == 1
            factory = factory[0]
            factoryName = factory.attrib['value']
            # rename factories to *_kokkos
            if factoryName in ('CoalesceDropFactory',
                               'CoarseMapFactory',
                               'CoordinatesTransferFactory',
                               'NullspaceFactory',
                               'SaPFactory',
                               'TentativePFactory',
                               'UncoupledAggregationFactory'):
                factory.attrib['value'] += '_kokkos'
            # delete FilteredAFactory
            if factoryName == 'FilteredAFactory':
                factoriesToDelete.append(factoryParams.attrib['name'])
                factoryList.remove(factoryParams)

        # delete references to deleted factories
        for factory in factoriesToDelete:
            for el in root.findall(".//*[@value='{}']/..".format(factory)):
                for el2 in el.findall(".//*[@value='{}']".format(factory)):
                    el.remove(el2)
    elif len(factoryList) == 0:
        # simple format
        useKokkos = ET.Element('Parameter')
        useKokkos.set('name', 'use kokkos refactor')
        useKokkos.set('type', 'bool')
        useKokkos.set('value', 'true')
        root.append(useKokkos)
    else:
        raise Exception("xml has wrong format.")

    tree.write(out)
