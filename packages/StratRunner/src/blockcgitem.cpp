/*
 * TODO
 * -figure out what the heck is up with the outputstream dealio
 *  	-It looks like the outputstream dealio has to do with where all the output from the solver goes (gee, who would have guessed that)
 *  	Maybe I just need to include the option to output the solver results to a file since that's the only other place I can think of to which a person would want to write output.
 * -can Timer label just be anything the user wants?
 *  	-I think it can be. I'm just going to add a new type of item called the FreeStringTreeItem to accomidate this.
 */

/*
 * blockcgitem.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
#include "blockcgitem.hpp"
#include "booltreeitem.hpp"
#include "stringtreeitem.hpp"
#include "inttreeitem.hpp"
#include "doubletreeitem.hpp"
#include "freestringtreeitem.hpp"

BlockCGItem::BlockCGItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, BlockCGItem::Type)
{
	setText(0,QString("BlockCG"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
}

void BlockCGItem::addParameters(){
	addChild(new BoolTreeItem("Adaptive Block Size", this, true));
	addChild(new BoolTreeItem("Show Maximum Residual Norm Only", this));

	addChild(new IntTreeItem("Block Size", "Size:", this, 1));
	addChild(new IntTreeItem("Maximum Iterations", "Maximum:", this, 1000));
	addChild(new IntTreeItem("Output Frequency", "Frequency:", this, -1));
	addChild(new IntTreeItem("Verbosity", "Verbosity:", this, 0));
	
	addChild(new DoubleTreeItem("Convergence Tolerance", "Tolerance:", this, 1e-08));
	addChild(new DoubleTreeItem("Orthogonalization Constant", "Constant:", this, -1));

	QStringList orthognalizationList = QStringList() << "DGKS" << "ICGS" << "IMGS";
	addChild(new StringTreeItem("Orthogonalization", "Orthogonalization:", orthognalizationList, this));

//	QStringList timerList = QStringList() << "Belos";
//	addChild(new StringTreeItem("Timer Label", "Label:", timerList, this));
	addChild(new FreeStringTreeItem("Timer Label", "Label:", this, "Belos"));
}

