/*
 * booltreeitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */
/*
 * TODO
 * -Change the way that the solver tree determines what the current value is. The optionOptions should really don't need to contain anything
 *
 */
#include "booltreeitem.hpp"
#include "parameterlisttreeitem.hpp"

BoolTreeItem::BoolTreeItem(QString name, ParameterTreeItem *parent, bool value)
	:ParameterTreeItem(parent, BoolTreeItem::Type)
{
	this->name = name;
	setText(0,name);
	setText(1,value ? "True" : "False");
	setText(2,QString("Bool"));
}

QStringList BoolTreeItem::getMenuOptions(){
	QStringList menuOptions;
	if(text(1) == "True"){
		menuOptions << "Set False";
	}
	else{
		menuOptions << "Set True";
	}
	return menuOptions;
}

QString BoolTreeItem::getOptionType(QString optionName){
	if(optionName == "Set False" || optionName == "Set True"){
		return "bool";
	}
	return "";
}
QString BoolTreeItem::getOptionLabel(QString optionName){
	return "";
}

QStringList BoolTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Set True"){
		optionOptions << "true";
	}
	else if(optionName == "Set False"){
		optionOptions << "false";
	}
	return optionOptions;
}

void BoolTreeItem::changeParameter(QString option, QString value){
	if(option == "Set True"){
		changeValue(true);
	}
	if(option == "Set False"){
		changeValue(false);
	}
}

void BoolTreeItem::changeValue(bool value){
	setText(1,value ? "True" : "False");
}
void BoolTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1) == "True" ? "1" : "0");
	xmlWriter.writeAttribute("type", "bool");
}

void BoolTreeItem::readInput(QXmlStreamReader &xmlReader){
	setText(1, xmlReader.attributes().value("value").toString() == "1" ? "True" : "False");
}

