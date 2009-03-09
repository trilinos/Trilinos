/*
 * belositem.cpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */
#include "belositem.hpp"
#include "belossolveritem.hpp"
#include "verboseobjectitem.hpp"

BelosItem::BelosItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, BelosSolverItem::Type)
{
	setText(0,"Belos");
	addParameterLists();
}

void BelosItem::addParameterLists(){
	addChild(new BelosSolverItem(this));
	addChild(new VerboseObjectItem(this));
}

