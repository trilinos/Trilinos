import xml.etree.ElementTree as et
import os

path = '.'

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

