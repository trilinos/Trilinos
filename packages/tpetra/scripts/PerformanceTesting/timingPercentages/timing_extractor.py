import xml.etree.ElementTree as et
import os
import getopt
import sys

def main():
        try:
                opts, args = getopt.getopt(sys.argv[1:], "p:", ["path="])
        except getopt.GetoptError as err:
                print(err)
                usage()
                sys.exit(2)
        path = '.'
        for k, v in opts:
                if k in ("-p", "--path"):
                        path = v
                else:
                        assert False, "unhandled option"
        addPercent(path)

def usage():
        print("please supply -p or --path arguments")


def addPercent(path) :
        for filename in os.listdir(path):
                if not filename.endswith('.xml'): continue
                fname = os.path.join(path, filename)
                xmlTree = et.parse(fname)
                first = True
                overallTime = 0.0
                for elem in xmlTree.iter():
                        if elem.tag == 'timing':
                                if first == True :
                                        overallTime = float(elem.attrib['value'])
                                        elem.set('percent', '100.0')
                                        first = False
                                else :			
                                        temp = float(elem.attrib['value'])
                                        pct = temp/overallTime * 100.0
                                        elem.set('percent', str(pct))

                xmlTree.write(fname)

if __name__ == "__main__":
        main()
