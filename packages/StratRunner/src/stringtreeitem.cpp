/*
 * stringtreeitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */
#include "stringtreeitem.hpp"
//#include "parameterlisttreeitem.hpp"

StringTreeItem::StringTreeItem(QString name, QString dialogLabel, QStringList validChoices, ParameterTreeItem *parent)
	:ParameterTreeItem(parent, StringTreeItem::Type)
{
	this->name = name;
	this->dialogLabel = dialogLabel;
	this->validChoices = validChoices;
	setText(0,name);
	setText(1,validChoices.at(0));
	setText(2,QString("String"));
}

QStringList StringTreeItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Value";
	return menuOptions;
}

QString StringTreeItem::getOptionType(QString optionName){
	if(optionName == "Change Value"){
		return "List";
	}
	return "";
}

QString StringTreeItem::getOptionLabel(QString optionName){
	return dialogLabel;
}

QStringList StringTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Change Value"){
		optionOptions = validChoices;
	}
	return optionOptions;
}

void StringTreeItem::changeParameter(QString option, QString value){
	if(option == "Change Value"){
		changeValue(value);
	}
}

void StringTreeItem::changeValue(QString value){
	setText(1,value);
}

void StringTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1));
	xmlWriter.writeAttribute("type", "string");
}

void StringTreeItem::readInput(QXmlStreamReader &xmlReader){
	if(validChoices.contains(xmlReader.attributes().value("value").toString())){
		setText(1, xmlReader.attributes().value("value").toString());
	}
}

