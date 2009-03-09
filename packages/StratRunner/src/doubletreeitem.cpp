/*
 * doubletreeitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */
/*
 * TODO
 * -Change the way that the solver tree determines what the current value is. The optionOptions should only contain the min and max value. The solver tree should just be able to look at the current value of the item using the text(1) function
 */
#include "doubletreeitem.hpp"
#include "parameterlisttreeitem.hpp"

DoubleTreeItem::DoubleTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, double value, QDoubleValidator *validator)
	:ParameterTreeItem(parent, DoubleTreeItem::Type)
{
	this->name = name;
	this->dialogLabel = dialogLabel;
	if(validator == 0){
			this->validator = new QDoubleValidator(treeWidget());
	}
	else{
		this->validator = validator;
	}
//	this->value = value;
	setText(0,name);
	setText(1,QString::number(value));
	setText(2,QString("Double"));
}

DoubleTreeItem::~DoubleTreeItem(){
//	delete validator;
}

QStringList DoubleTreeItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Value";
	return menuOptions;
}

QString DoubleTreeItem::getOptionType(QString optionName){
	if(optionName == "Change Value"){
		return "double";
	}
	return "";
}

QString DoubleTreeItem::getOptionLabel(QString optionName){
	return dialogLabel;
}

QStringList DoubleTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Change Value"){
//		optionOptions << QString::number(value) <<  QString::number(validator->bottom()) << QString::number(validator->top());
		optionOptions << text(1) <<  QString::number(validator->bottom()) << QString::number(validator->top());
	}
	return optionOptions;
}

void DoubleTreeItem::changeParameter(QString option, QString value){
	if(option == "Change Value"){
		changeValue(value.toDouble());
	}
}

void DoubleTreeItem::changeValue(double value){
	//this->value = value;
	setText(1,QString::number(value));
}
/*
double DoubleTreeItem::getValue(){
		//return value;
	return text(1).toDouble();
}

double DoubleTreeItem::getMin(){
	return validator->bottom();
}

double DoubleTreeItem::getMax(){
	return validator->top();
}*/

void DoubleTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1));
	xmlWriter.writeAttribute("type", "double");
}

void DoubleTreeItem::readInput(QXmlStreamReader &xmlReader){
	setText(1, xmlReader.attributes().value("value").toString());
}

