from optparse import OptionParser

def remove_method_body(arg, level):
    chars = ""
    n = 0
    ignoreNext = False

    for c in arg:
        if ignoreNext:
            ignoreNext = False
            if n < level:
                chars += c
        else:
            if c == '@':
                ignoreNext = True

            if c == '{':
                n += 1
                if n == level:
                    chars += ";"

            if n < level:
               chars += c

            if c == '}':
               n -= 1

    return chars

parser = OptionParser()
parser.add_option("-i", "--input",  dest="i", help="input file",  metavar="FILE")
parser.add_option("-o", "--output", dest="o", help="output file", metavar="FILE")
parser.add_option("-l", "--level",  dest="l", help="level", type="int")

(options, args) = parser.parse_args()

str = open(options.i, 'r').read()
o = remove_method_body(str, options.l)

open(options.o, 'w').write(o)
