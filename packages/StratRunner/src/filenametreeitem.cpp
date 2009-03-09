/*
 * filenametreeitem.cpp
 *
 *  Created on: Dec 19, 2008
 *      Author: kurtis
 */
#include "filenametreeitem.hpp"

FileNameTreeItem::FileNameTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent)
	:ParameterTreeItem(parent, FileNameTreeItem::Type)
{
	this->name = name;
	this->dialogLabel = dialogLabel;
//	this->value = "";
	setText(0,name);
//	setText(1,value);
	setText(1,"none");
	setText(2,QString("String"));
}

QStringList FileNameTreeItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Value";
	return menuOptions;
}

QString FileNameTreeItem::getOptionType(QString optionName){
	if(optionName == "Change Value"){
		return "filename";
	}
	return "";
}

QString FileNameTreeItem::getOptionLabel(QString optionName){
	return dialogLabel;
}

QStringList FileNameTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Change Value"){
		//optionOptions.append(value);
		optionOptions.append(text(1));
	}
	return optionOptions;
}

void FileNameTreeItem::changeParameter(QString option, QString value){
	if(option == "Change Value"){
		changeValue(value);
	}
}

void FileNameTreeItem::changeValue(QString value){
	//this->value = value;
	setText(1,value);
}

void FileNameTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1));
	xmlWriter.writeAttribute("type", "string");
}

void FileNameTreeItem::readInput(QXmlStreamReader &xmlReader){
	setText(1, xmlReader.attributes().value("value").toString());
}

