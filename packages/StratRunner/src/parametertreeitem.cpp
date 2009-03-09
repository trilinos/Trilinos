/*
 * parametertreeitem.cpp
 *
 *  Created on: Dec 6, 2008
 *      Author: kurtis
 */
#include "parametertreeitem.hpp"
#include "solvertree.hpp"

ParameterTreeItem::ParameterTreeItem(SolverTree *parent, int type)
	:QTreeWidgetItem(parent, type)
{

}

ParameterTreeItem::ParameterTreeItem(ParameterTreeItem *parent, int type)
	:QTreeWidgetItem(parent, type)
{

}

ParameterTreeItem* ParameterTreeItem::findChildByName(QString name){
	for(int i=0; i<childCount(); i++){
		if(child(i)->text(0) == name){
			return dynamic_cast<ParameterTreeItem*>(child(i));
		}
	}
	return 0;
}

