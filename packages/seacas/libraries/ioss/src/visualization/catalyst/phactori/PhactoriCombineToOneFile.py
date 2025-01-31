
#combine phactori.py and all operations (and later other pieces) into
#phactori_onefile.py, and then combine phactori_onefile.py and
#PhactoriDriver_modular.py into PhactoriDriver_onefile.py

def ReplaceQuotesWithEscapes(inStr):

    outStr = ""
    for ndx in range(0,len(inStr)):
        ch = inStr[ndx]
        if ch == "\n":
            outStr = outStr + "\\n"
        elif ch == '"':
            outStr = outStr + "'"
        else:
            outStr = outStr + ch

    return outStr

def RemoveTripleQuotes(inlines):
    outlines = []
    ndx = 0
    isDone = False
    while not isDone:
        if ndx >= len(inlines):
            isDone = True
            continue

        #copy lines until you see """
        if inlines[ndx].find('"""') == -1:
            outlines.append(inlines[ndx])
            ndx += 1
            continue

        oneLine = inlines[ndx]
        qndx1 = oneLine.find('"""')
        qndx2 = oneLine[qndx1+3:].find('"""')

        #if this line is a comment, just pass it through
        commentndx = inlines[ndx].find('#')
        if commentndx != -1:
            if commentndx < qndx1:
                outlines.append(inlines[ndx])
                ndx += 1
                continue

        print(("found tq:", ndx))
        print(oneLine)

        print((qndx1, qndx2))
        #handle case where line contains two """
        if qndx2 > -1:
            insideQuotes = oneLine[qndx1+3:qndx1+3+qndx2]
            insideQuotesWithEscapes = ReplaceQuotesWithEscapes(insideQuotes)
            chgLine = oneLine[0:qndx1] + '"' + insideQuotesWithEscapes + '"' + oneLine[qndx1+3+qndx2+3:]
            outlines.append(chgLine)
            ndx += 1
            print(("fixedone1", ndx))
            print("complete replacement line 2:")
            print(chgLine)
            continue

        print("found2")
        print(oneLine)
        print((qndx1, qndx2))

        prefix = oneLine[0:qndx1]

        chgLine = oneLine[qndx1+3:]

        #find next line that has """, closing quote
        qend_ndx = ndx + 1
        while inlines[qend_ndx].find('"""') == -1:
            chgLine = chgLine + inlines[qend_ndx]
            qend_ndx += 1

        oneLine = inlines[qend_ndx]

        print(("end of triple quote", qend_ndx))
        print(oneLine)

        #freak out if this line has two """
        qndx1 = oneLine.find('"""')
        qndx2 = oneLine[qndx1+3:].find('"""')
        if qndx2 != -1:
            print("can't handle close triple quote and opening triple quote on same line\n")
            exit(-100)

        chgLine = chgLine + oneLine[:qndx1]
        postfix = oneLine[qndx1+3:]

        chgLine = ReplaceQuotesWithEscapes(chgLine)
        chgLine = prefix + '"' + chgLine + '"' + postfix

        print("complete replacement line:")
        print(chgLine)
        outlines.append(chgLine)
        ndx = qend_ndx + 1
    return outlines



def GetCombineSectionFromClass(inClassName):
    """get lines from subclass (including subdirectory name)and return
       to be added to merged file"""
    fname = inClassName + ".py"
    pp = open(fname, "r")
    if pp == None:
        print(("could not open file:", fname))
        exit(-5)
    pplines = pp.readlines()
    pp.close()
    #find beginning/end indexes
    ndx = 0
    beginNdx = -1
    endNdx = -1
    while ndx < len(pplines):
        if pplines[ndx].startswith("#phactori_combine_to_single_python_file_subpiece_begin_1"):
            beginNdx = ndx
        if pplines[ndx].startswith("#phactori_combine_to_single_python_file_subpiece_end_1"):
            endNdx = ndx
            break
        ndx += 1

    if (beginNdx < 0) or (endNdx < 0):
        print(("file:", fname))
        if beginNdx < 0:
            print("could not find #phactori_combine_to_single_python_file_subpiece_begin_1")
        if endNdx < 0:
            print("could not find #phactori_combine_to_single_python_file_subpiece_end_1")
        exit(-6)

    retList = pplines[beginNdx:endNdx+1]
    return retList


ff = open("phactori.py","r")
phlines = ff.readlines()
oo = open("phactori_onefile.py", "w")

isDone = False
phndx = 0
while(isDone == False):
    #copy until #phactori_combine_to_single_python_file_parent_ line
    while (phlines[phndx].startswith("#phactori_combine_to_single_python_file_parent_") == False):
        oo.write(phlines[phndx])
        phndx += 1
        if phndx >= len(phlines):
            isDone = True
            break

    if isDone == True:
        break

    #go ahead and copy combine comment
    oo.write(phlines[phndx])
    phndx += 1

    #make sure next line starts with from classname import *
    if (phlines[phndx].startswith("from ") == False):
        print(("expecting 'from X import *' line to replace, ndx", phndx))
        exit(-1)

    #get class to combine in
    xx = phlines[phndx].split(" ")

    if xx[2] != "import":
        print(("expecting (2) 'from X import *' line to replace, ndx", phndx))
        exit(-2)

    if xx[3][0] != "*":
        print(("expecting (3) 'from X import *' line to replace, ndx", phndx))
        exit(-3)

    zz = xx[1].split(".")
    if len(zz) != 2:
        print(("expecting something like 'from import Operation.PhactoriGroupOperation import*', ndx", phndx))
        exit(-4)

    subdir = zz[0]
    combineFileClass = zz[1]

    print(("combining dir, class, ndx:", subdir, combineFileClass, phndx))

    #add in import line with comment
    cmtimprtln = "#" + phlines[phndx]
    oo.write(cmtimprtln)
    phndx += 1

    #add in correct segment of external class (get correct lines from
    #sub piece and merge them in
    targetFilename = subdir + "/" + combineFileClass
    classToCombineAsStringList = GetCombineSectionFromClass(targetFilename)
    for ssl in classToCombineAsStringList:
      oo.write(ssl)

oo.close()

#now combine into one PhactoriDriver file

pof = open("phactori_onefile.py", "r")
poflineswtq = pof.readlines()
pof.close()
poflines = RemoveTripleQuotes(poflineswtq)

phdrv = open("PhactoriDriver_modular.py", "r")
if phdrv == None:
    print(("could not open file:", "PhactoriDriver_modular.py"))
    exit(-8)
phdrvlines = phdrv.readlines()
phdrv.close()


phdrvof = open("PhactoriDriver_onefile.py", "w")
if phdrvof == None:
    print(("could not open file:", "PhactoriDriver_onefile.py"))
    exit(-9)

isDone = False
lnndx = 0
frmphimprtlnfound = False
while(isDone == False):
    #copy until "from phactori import *" line
    while (phdrvlines[lnndx].startswith("from phactori import *") == False):
        phdrvof.write(phdrvlines[lnndx])
        lnndx += 1
        if lnndx >= len(phdrvlines):
            isDone = True
            break

    if isDone == True:
        break

    if phdrvlines[lnndx].startswith("from phactori import *"):
        if frmphimprtlnfound:
            print("error: from phactori import * in in twice!\n")
            exit(-11)
        frmphimprtlnfound = True

    phdrvof.write("#" + phdrvlines[lnndx])
    phdrvof.write("######phactori.py inserted here by PhactoriCombineToOneFile.py\n\n")
    lnndx += 1
    for oneLine in poflines:
        phdrvof.write(oneLine)

phdrvof.close()

if frmphimprtlnfound == False:
    print("error, from phactori import * not found in PhactoriDriver_modular.py")
    exit(-12)
