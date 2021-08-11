<!DOCTYPEhtml>
 <html lang="en-US">
  <body>

<h1>Factory Graph</h1>
<p>A python tool to combine files into a nested dot file. All factories that
are in each file will be placed into a graph. That graph will be saved to a
.dot file and a .png file.</p>

<h2>Rules</h2>
<ol>
	<li>Some nodes will appear on multiple levels. Those nodes should only 
	appear on the first (lowest number) subgraph in which they exist (and the
	arrows should route between subgraphs).</li>
	<li>"NoFactory" is a special case node and exists on every subgraph. Unlike other duplicates, that should just remain on each individual subgraph.</li>
	<li>Nodes with no edges can be removed.</li>
</ol>

<h2>How to Use</h2>
<ol>
	<li>The python file <b>MUST</b> take at least one comman line argument. Each argument <b>MUST</b> be a .dot file. For example, in the terminal type: <code>python3 mergeLevelGraphs.py file.dot file1.dot</code>.</li>
	<li>There is a optional command line argument <code>-o</code>. If used the user <b>MUST</b> pass two files paths that the outputted .dot file and .png file will be saved to, for example, <code>python3 mergeLevelGraphs.py file.dot file1.dot -o fileLocation2SaveDot.dot fileLocation2SavePng.png</code>. If the optional -o command is not used then mergeLevelGraphs.py will write to output.dot and output.png. If those files exist then they will be overwritten.</li>
	<li>Create an instance of the class, for example, <code>merge = MergeLevelGraphs()</code>.</li>
	<li>To merge all the factories call the merge method, for example, <code>merge.merge()</code>.</li>
	<li>To see either the subgraph or master graph use the properties <code>merge.getSubGraph</code> and <code>merge.getMasterGraph</code>.</li>
	<li>To save to a dot file and png call the write method, for example, <code>merge.write()</code>.
</ol>

 </html>
</body>
