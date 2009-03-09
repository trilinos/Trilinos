/*
 * preconitem.cpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */

#include "preconitem.hpp"
#include "nopreconitem.hpp"
#include "mlitem.hpp"
#include "ifpackitem.hpp"

PreconItem::PreconItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, PreconItem::Type)
{
	setText(0,"Preconditioner Type");
	setText(1,"None");
	setText(2,"String");

	preconChild = new NoPreconItem(this);

	preconsList << "ML" << "Ifpack";
}

QStringList PreconItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Preconditioner";
	return menuOptions;
}

QString PreconItem::getOptionType(QString optionName){
	return "List";
}

QString PreconItem::getOptionLabel(QString optionName){
	return "Preconditioner:";
}

QStringList PreconItem::getOptionOptions(QString optionName){
	return preconsList;
}

void PreconItem::changeParameter(QString optionName, QString value){
	if(value != text(1)){
		preconsList << text(1);
		if(value == "None"){
			delete preconChild;
			preconChild = new NoPreconItem(this);
			addChild(preconChild);
		}
		else if(value == "ML"){
			delete preconChild;
			preconChild = new MLItem(this);
			addChild(preconChild);
		}
		else if(value == "Ifpack"){
			delete preconChild;
			preconChild = new IfpackItem(this);
			addChild(preconChild);
		}
		preconsList.removeAll(value);
		setText(1,value);
	}
}

void PreconItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name",text(0));
	xmlWriter.writeAttribute("type",text(2));
	xmlWriter.writeAttribute("value",text(1));
	xmlWriter.writeStartElement("ParameterList");
	xmlWriter.writeAttribute("name", "Preconditioner Types");
	preconChild->writeOutput(xmlWriter);
	xmlWriter.writeEndElement();
}

void PreconItem::readInput(QXmlStreamReader &xmlReader){
	changeParameter("Change Solver", xmlReader.attributes().value("value").toString());	
	ParameterListTreeItem::readInput(xmlReader);
}

