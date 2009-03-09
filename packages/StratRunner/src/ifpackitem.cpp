/*
 * ifpackitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

#include "ifpackitem.hpp"
#include "inttreeitem.hpp"
#include "stringtreeitem.hpp"
#include "ifpacksettingsitem.hpp"
#include "verboseobjectitem.hpp"


IfpackItem::IfpackItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, IfpackItem::Type)
{
	setText(0,QString("Ifpack"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
}

void IfpackItem::addParameters(){
	addChild(new IntTreeItem("Overlap", "Overlap:", this));

	QStringList preconditionTypes = QStringList() << "point relaxation" << "point relaxation stand-alone" << "block relaxation" << "block relaxation stand-alone" << "block relaxation stand-alone (ILU)" << "block relaxation stand-alone (Amesos)" << "block relaxation (Amesos)" << "Amesos" << "Amesos stand-alone" << "IC" << "IC stand-alone" << "ICT" << "ICT stand-alone" << "ILU" << "ILU stand-alone" << "ILUT" << "ILUT stand-alone" << "Chebyshev";
	addChild(new StringTreeItem("Prec Type", "Type:", preconditionTypes, this));
}

void IfpackItem::addParameterLists(){
	addChild(new IfpackSettingsItem(this));
	addChild(new VerboseObjectItem(this));
}

