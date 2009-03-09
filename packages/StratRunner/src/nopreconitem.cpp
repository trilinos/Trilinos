/*
 * nopreconitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */
#include "nopreconitem.hpp"

NoPreconItem::NoPreconItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, NoPreconItem::Type)
{
	setText(0,QString("No Preconditioner"));
	setText(1,QString(""));
	setText(2,QString(""));
}

