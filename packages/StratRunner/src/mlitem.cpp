/*
 * mlitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */
#include "mlitem.hpp"
#include "stringtreeitem.hpp"
#include "mlsettingsitem.hpp"

MLItem::MLItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, MLItem::Type)
{
	setText(0,QString("ML"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void MLItem::addParameters(){
	QStringList baseMethodDefaults = QStringList() << "none" << "SA" << "DD" << "DD-ML" << "maxwell";
	addChild(new StringTreeItem("Base Method Defaults", "Defaults:", baseMethodDefaults, this));
}

void MLItem::addParameterLists(){
	addChild(new MLSettingsItem(this));
}

