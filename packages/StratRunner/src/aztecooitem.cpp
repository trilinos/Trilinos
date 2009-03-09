/*
 * aztecooitem.cpp
 *
 *  Created on: Dec 5, 2008
 *      Author: kurtis
 */
#include "aztecooitem.hpp"
#include "booltreeitem.hpp"
#include "inttreeitem.hpp"
#include "doubletreeitem.hpp"
#include "aztecoosettingsitem.hpp"
#include "verboseobjectitem.hpp"

AztecOOItem::AztecOOItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, AztecOOItem::Type)
{
	setText(0,QString("AztecOO"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void AztecOOItem::addParameters(){
	addChild(new BoolTreeItem("Output Every RHS", this));
}

void AztecOOItem::addParameterLists(){
	QList<ParameterTreeItem*> adjointItems;
	adjointItems.append(new IntTreeItem("Max Iterations", "Max:", this, 400));
	adjointItems.append(new DoubleTreeItem("Tolerance", "Tolerance:", this, 1e-06));	
	adjointItems.append(new AztecOOSettingsItem(this));
	addChild(new ParameterListTreeItem("Adjoint Solve", this, adjointItems));

	QList<ParameterTreeItem*> forwardItems;
	forwardItems.append(new IntTreeItem("Max Iterations", "Max:", this, 400));
	forwardItems.append(new DoubleTreeItem("Tolerance", "Tolerance:", this, 1e-06));	
	forwardItems.append(new AztecOOSettingsItem(this));
	addChild(new ParameterListTreeItem("Forward Solve", this, forwardItems));
	addChild(new VerboseObjectItem(this)); 
}

