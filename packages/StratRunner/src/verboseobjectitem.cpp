/*
 * verboseobjectitem.cpp
 *
 *  Created on: Dec 19, 2008
 *      Author: kurtis
 */
#include "verboseobjectitem.hpp"
#include "stringtreeitem.hpp"
#include "filenametreeitem.hpp"

VerboseObjectItem::VerboseObjectItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, VerboseObjectItem::Type)
{
	setText(0,QString("Verbose Object"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void VerboseObjectItem::addParameters(){
	addChild(new FileNameTreeItem("Output File", "File Name:", this));
	QStringList verbosityChoices = QStringList() << "default" << "none" << "low" << "medium" << "high" << "extreme";
	addChild(new StringTreeItem("Verbosity Level", "Level", verbosityChoices, this));
}


