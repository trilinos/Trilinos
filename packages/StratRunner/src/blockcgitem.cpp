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

	addChild(new FreeStringTreeItem("Timer Label", "Label:", this, "Belos"));
}

