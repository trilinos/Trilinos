import os
import sys
import pydot

class MergeLevelGraphs:
	'''
	A python tool to combine files into a nested dot file. All factories that
	are in each file will be placed into a graph. That graph will be saved to a
	.dot file and a .png file.\n
	Rules:
	1) Some nodes will appear on multiple levels. Those nodes should only 
	appear on the first (lowest number) subgraph in which they exist (and the
	arrows should route between subgraphs).\n
	2) "NoFactory" is a special case node and exists on every subgraph. Unlike
	other duplicates, that should just remain on each individual subgraph.\n
	3) Nodes with no edges can be removed.
	'''
	__graphs = {}
	# store each subgraph in the master graph
	__masterGraph = pydot.Dot(graph_name="Master Graph",
		suppress_disconnected=True, simplify=True)
	__nodeMap = {}
	__repeatNodeMap = {}
	__graphCnt = 0

	def __init__(self):
		'''Convert each .dot file to a graph and append to graph list'''
		if len(sys.argv) < 1:
			raise ValueError("0 command line argument passed, but expected at least 1")

		graphLst = []

		if "-o" in sys.argv:
			graphLst = sys.argv[1:-3]
		else:
			graphLst = sys.argv[1:]

		for ii in graphLst:
			isContinue = False

			while os.path.exists(ii) == False:
				ans = input(f"File {ii} does not exists. Enter \"y\" to replace it: ")

				if ans.lower().strip() == "y":
					ii = input("Enter file location: ")
				else:
					isContinue = True
					break
			
			if isContinue:
				continue

			index = ""
			if sys.platform == "win32":
				index = ii.rindex("\\") + 1
			else:
				index = ii.rindex("/") + 1

			fileName = ii[index: ii.rindex(".")]

			if fileName in self.__graphs.values():
				raise ValueError("File names cannot be the same, even if path is different:", fileName)

			self.__graphs[pydot.graph_from_dot_file(ii)[0]] = fileName

	def __deleteNodes(self) -> None:
		'''Delete nodes that are no longer needed.'''
		for point in self.__delNodes:
			if ((str(self.__graphCnt) + "_") in point) == False:
				continue

			self.__subG.del_node(point)

			for edge in self.__subG.get_edge_list():
				POINTS = edge.obj_dict["points"]

				if point in POINTS:
					self.__subG.del_node(POINTS)

	def __findNode(self, NAME: str) -> list:
		'''
		Find the label that is not in the subgraph by searching masterGraph for it.

		# Params:
		NAME - name of the node being searched for
		# Returns:
		A list that contains the node that was being searched for
		'''
		num = "0"

		if "_" in NAME: # get the correct cluster number
			num = NAME[:NAME.index("_")].strip("\"").strip("'").strip()

		# DON'T use obj_dict becuase it won't return a list with the node object
		return self.__masterGraph.get_subgraph(f"cluster_{num}")[0].get_node(NAME)

	def __addEdge(self, SRC_DST: tuple, ATTRIBUTES: tuple) -> None:
		'''
		Add edges to the subgraph.

		# Params:
		SRC_DST - index 0 is SRC and 1 is DST\n
		ATTRIBUTES - 0 is color 1 is label
		'''
		self.__subG.add_edge(pydot.Edge(SRC_DST[0], SRC_DST[1],
			color=ATTRIBUTES[0], label=ATTRIBUTES[1]))

	def __replaceConnection(self, LABEL: str, ATTRIB_DICT: dict, IS_0_SRC: bool=True) -> None:
		'''
		Delete node edges that need to be replaced. Then, determin if
		self.__NODE_0 is source and replace deleted edge the node.

		# Params:
		LABEL - Label of the node that be used for replacement\n
		ATTRIB_DICT - key: type of attribute, value: attribute\n
		IS_0_SRC - If self.__NODE_0 is src replace that node else replace
		self.__NODE_1
		'''
		self.__subG.del_edge(self.__NODE_0[0], self.__NODE_1[0])

		if IS_0_SRC:
			self.__addEdge((self.__repeatNodeMap[LABEL][0], self.__NODE_1[0]),
				(ATTRIB_DICT["color"], ATTRIB_DICT["label"]))
		else:
			self.__addEdge((self.__NODE_0[0], self.__repeatNodeMap[LABEL][0]),
				(ATTRIB_DICT["color"], ATTRIB_DICT["label"]))
	
	def __renameNodes(self):
		'''Rename all the nodes and remap their subgraphs.'''
		KEYS = list(self.__subG.obj_dict["nodes"])

		for ii in range(len(KEYS)):
			NODE_DATA = self.__subG.obj_dict["nodes"][KEYS[ii]][0]
			# delete old node
			self.__subG.del_node(KEYS[ii])
			# create new node to replace old node
			newName = f"{self.__graphCnt}_{KEYS[ii]}" # fileNumber_nodeName
			self.__subG.add_node(pydot.Node(name=newName,
				label=NODE_DATA["attributes"]["label"]))

			# remap the node edges
			for edge in self.__subG.get_edge_list():
				POINTS = edge.obj_dict["points"]
				
				# len(POINTS) WILL ALWAYS be 2, so Big O won't be that bad
				for point in POINTS:
					if KEYS[ii] != point:
						continue

					ATTRIB_DICT = edge.obj_dict["attributes"]
					ATTRIB = (ATTRIB_DICT["color"], ATTRIB_DICT["label"])

					# replace the old node's connections
					if KEYS[ii] == POINTS[0]:
						self.__addEdge((newName, POINTS[1]),
							(ATTRIB[0], ATTRIB[1]))
					else:
						self.__addEdge((POINTS[0], newName),
							(ATTRIB[0], ATTRIB[1]))
					self.__subG.del_edge(POINTS)

	def __renameSubgraphs(self) -> None:
		'''Change each subgraph name from cluster_# to cluster_fileName.'''
		cnt = 0

		for value in self.__graphs.values():
			NAME = f"cluster_{cnt}"
			NEW_NAME = f"cluster_{value}"

			self.__masterGraph.obj_dict["subgraphs"][NEW_NAME] = [self.__masterGraph.obj_dict["subgraphs"][NAME].pop()]
			self.__masterGraph.obj_dict["subgraphs"][NEW_NAME][0]["name"] = NEW_NAME
			del self.__masterGraph.obj_dict["subgraphs"][NAME]
			cnt += 1

	def merge(self) -> None:
		'''Put all factories in order according to the rules'''
		for graph in self.__graphs:
			graph.obj_dict["type"] = "subgraph"
			graph.obj_dict["name"] = f"cluster_{self.__graphCnt}"
			graph.set_suppress_disconnected(True)
			graph.set_label(f"cluster_{self.__graphs[graph]}")

			self.__subG = pydot.Subgraph(graph_name=f"cluster_{self.__graphCnt}",
				obj_dict=graph.obj_dict, suppress_disconnected=True)

			if self.__graphCnt == 0:
				for key in self.__subG.obj_dict["nodes"]:
					self.__nodeMap[self.__subG.obj_dict["nodes"][key][0]["attributes"]
						["label"]] = self.__subG.get_node(key)

				self.__graphCnt += 1
				self.__masterGraph.add_subgraph(self.__subG)
				continue

			self.__renameNodes()

			self.__delNodes = []
			for key in self.__subG.obj_dict["nodes"]:
				LABEL = self.__subG.obj_dict["nodes"][key][0]["attributes"]["label"]

				# if node in subgraph is a repeat node, map it and remove it from subG
				if(((LABEL in self.__nodeMap) == False) or ("NoFactory" in LABEL)):
					self.__nodeMap[LABEL] = self.__subG.get_node(key)
					continue

				edgeLst = self.__subG.get_edge_list()
				searchPoint = ""

				if (LABEL in self.__repeatNodeMap) == False:
					self.__repeatNodeMap[LABEL] = self.__nodeMap[LABEL]

				for edge in edgeLst:
					# map the nodes to their new connections
					POINTS = edge.obj_dict["points"]

					if((searchPoint != "") and ((searchPoint in POINTS) == False)):
						continue

					self.__NODE_0 = self.__subG.get_node(POINTS[0])
					self.__NODE_1 = self.__subG.get_node(POINTS[1])

					if len(self.__NODE_0) < 1:
						self.__NODE_0 = self.__findNode(POINTS[0])
					elif len(self.__NODE_1) < 1:
						self.__NODE_1 = self.__findNode(POINTS[1])

					# connect repeat nodes to it's counter part
					if LABEL == self.__NODE_0[0].obj_dict["attributes"]["label"]:
						searchPoint = POINTS[0]
						self.__replaceConnection(LABEL,
							edge.obj_dict["attributes"])

					elif LABEL == self.__NODE_1[0].obj_dict["attributes"]["label"]:
						searchPoint = POINTS[1]
						self.__replaceConnection(LABEL,
							edge.obj_dict["attributes"], False)

					if searchPoint != "":
						if (searchPoint in self.__delNodes) == False: 
							self.__delNodes.append(searchPoint)

			self.__deleteNodes()
			# add subgraph to masterGraph
			self.__masterGraph.add_subgraph(self.__subG)
			self.__graphCnt += 1
		self.__renameSubgraphs()

	def __getFileName(self, PATH: str, ext: str=".dot") -> int:
		'''
		Finds a file number that does not exists and sets fileCnt to it.

		# Params:
		PATH - file path to save the .dot/.png file\n
		fileName - the file extension determin if an image file
		exists or a dot file exists\n
		# Return:
		file number that does not exist
		'''
		fileCnt = 0

		while os.path.exists(os.path.join(PATH, f"{fileCnt}{ext}")):
			fileCnt += 1
		return fileCnt

	def write(self) -> None:
		'''Save each subgraph to a file.'''
		PATH = os.path.dirname(__file__)
		graphFile = ""
		imageFile = ""

		if "-o" in sys.argv:
			try:
				graphFile = sys.argv[-2]
			except:
				raise RuntimeError("Pass two file paths in the command line arguments to save to file")

			imageFile = sys.argv[-1]

			if os.access(graphFile, os.X_OK) == False:
				if graphFile[-4:] != ".dot":
					raise RuntimeError(f"File must be a .dot file {graphFile}")

			if os.access(imageFile, os.X_OK) == False:
				if imageFile[-4:] != ".png":
					raise RuntimeError(f"File must be a .png file {imageFile}")

			while os.path.exists(graphFile):
				ans = input(f"Are you sure you want to overwrite this file: {graphFile}\n Enter \"y\" to pick a different file name: ")
				if ans.strip().lower() == "y":
					ans = input("Enter file path: ")
					graphFile = ans
				else:
					break

			while os.path.exists(imageFile):
				ans = input(f"Are you sure you want to overwrite this file: {imageFile}\n Enter \"y\" to pick a different file name: ")
				if ans.strip().lower() == "y":
					ans = input("Enter file path: ")
					imageFile = ans
				else:
					break
		else:
			graphFile = "output.dot"
			imageFile = "output.png"

		try:
			self.__masterGraph.write(graphFile)
		except:
			raise(f"There is an issue saving file path: {graphFile}")

		try:
			self.__masterGraph.write_png(imageFile)
		except:
			raise(f"There is an issue saving to file path: {imageFile}")

	def __getMasterGraph(self) -> pydot.Dot:
		'''
		# Returns:
		An instance of pydot.Dot that contains all merged
		subgraphs
		'''
		return self.__masterGraph

	def __getSubGraph(self) -> pydot.Subgraph:
		'''
		# Returns:
		An instance of pydot.Subgraph that contains
		all nodes, edges, and etc in the subgraph. Depending on when this
		method is called the subgraph may or may not be 100% ordered yet. It is
		recommended to call MergeLevelGraphs.getMasterGraph() instead and use
		MergeLevelGraphs.getMasterGraph.get_subgraph(name) to find the
		appropriate subgraph.
		'''
		return self.__subG

	getMasterGraph = property(__getMasterGraph)
	getSubGraph = property(__getSubGraph)

if __name__ == "__main__":
	merge = MergeLevelGraphs()
	merge.merge()
	merge.write()
