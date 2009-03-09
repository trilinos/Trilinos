/*
 * freestringtreeitem.cpp
 *
 *  Created on: Feb 10, 2009
 *      Author: kurtis
 */

#include "freestringtreeitem.hpp"
#include "parameterlisttreeitem.hpp"

FreeStringTreeItem::FreeStringTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, QString defaultValue)
	:ParameterTreeItem(parent, FreeStringTreeItem::Type)
{
	this->name = name;
	this->dialogLabel = dialogLabel;
	setText(0,name);
	setText(1,defaultValue);
	setText(2,QString("String"));
}

QStringList FreeStringTreeItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Value";
	return menuOptions;
}

QString FreeStringTreeItem::getOptionType(QString optionName){
	if(optionName == "Change Value"){
		return "freestring";
	}
	return "";
}

QString FreeStringTreeItem::getOptionLabel(QString optionName){
	return dialogLabel;
}

QStringList FreeStringTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Change Value"){
		optionOptions << "";
	}
	return optionOptions;
}

void FreeStringTreeItem::changeParameter(QString option, QString value){
	if(option == "Change Value"){
		changeValue(value);
	}
}

void FreeStringTreeItem::changeValue(QString value){
	setText(1,value);
}

void FreeStringTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1));
	xmlWriter.writeAttribute("type", "string");
}

void FreeStringTreeItem::readInput(QXmlStreamReader &xmlReader){
	setText(1, xmlReader.attributes().value("value").toString());
}

