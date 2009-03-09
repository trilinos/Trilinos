/*
 * TODO
 * -figure out what the heck is up with the outputstream dealio
 *  	-see blockcg for solution
 * -can Timer label just be anything the user wants?
 *  	-see blockcg for solution
 */

/*
 * psuedoblockgmresitem.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
#include "psuedoblockgmresitem.hpp"
#include "booltreeitem.hpp"
#include "inttreeitem.hpp"
#include "stringtreeitem.hpp"
#include "doubletreeitem.hpp"
#include "freestringtreeitem.hpp"

PsuedoBlockGMRESItem::PsuedoBlockGMRESItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, PsuedoBlockGMRESItem::Type)
{
	setText(0,QString("PsuedoBlockGMRES"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
}

void PsuedoBlockGMRESItem::addParameters(){
	addChild(new BoolTreeItem("Adaptive Block Size", this, true));
	addChild(new BoolTreeItem("Show Maximum Residual Norm Only", this));

	addChild(new IntTreeItem("Block Size", "Size:", this, 1));
	addChild(new IntTreeItem("Maximum Iterations", "Maximum:", this, 1000));
	addChild(new IntTreeItem("Maximum Restarts", "Maximum:", this, 20));
	addChild(new IntTreeItem("Output Frequency", "Frequency:", this, -1));
	addChild(new IntTreeItem("Verbosity", "Verbosity:", this, 0));
	
	addChild(new DoubleTreeItem("Convergence Tolerance", "Tolerance:", this, 1e-08));
	addChild(new DoubleTreeItem("Orthogonalization Constant", "Constant:", this, -1));

	QStringList orthognalizationList = QStringList() << "DGKS" << "ICGS" << "IMGS";
	addChild(new StringTreeItem("Orthogonalization", "Orthogonalization:", orthognalizationList, this));

	//QStringList explicitScailingList = QStringList() << "
	//QStringList implicitScailingList = QStringList() <<	
	//QStringList timerList = QStringList() << "Belos";
	//addChild(new StringTreeItem("Timer Label", "Label:", timerList, this));
	addChild(new FreeStringTreeItem("Timer Label", "Label:", this, "Belos"));
}

