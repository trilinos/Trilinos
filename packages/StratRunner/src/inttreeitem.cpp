/*
 * inttreeitem.cpp
 *
 *  Created on: Dec 8, 2008
 *      Author: kurtis
 */

/*
 * TODO
 * -Change the way that the solver tree determines what the current value is. The optionOptions should only contain the min and max value. The solver tree should just be able to look at the current value of the item using the text(1) function
 */
#include "inttreeitem.hpp"
#include "parameterlisttreeitem.hpp"

IntTreeItem::IntTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent, int value, QIntValidator *validator)
	:ParameterTreeItem(parent, IntTreeItem::Type)
{
	this->name = name;
	this->dialogLabel = dialogLabel;
	if(validator == 0){
			this->validator = new QIntValidator(treeWidget());
	}
	else{
		this->validator = validator;
	}
	setText(0,name);
	setText(1,QString::number(value));
	setText(2,QString("Int"));
}

IntTreeItem::~IntTreeItem(){
//	delete validator;
}

QStringList IntTreeItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Value";
	return menuOptions;
}

QString IntTreeItem::getOptionType(QString optionName){
	if(optionName == "Change Value"){
		return "int";
	}
	return "";
}

QString IntTreeItem::getOptionLabel(QString optionName){
	return dialogLabel;
}

QStringList IntTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	if(optionName == "Change Value"){
		optionOptions << text(1)  << QString::number(validator->bottom()) << QString::number(validator->top());
	}
	return optionOptions;
}

void IntTreeItem::changeParameter(QString option, QString value){
	if(option == "Change Value"){
		changeValue(value.toInt());
	}
}

void IntTreeItem::changeValue(int value){
	setText(1,QString::number(value));
}
/*
int IntTreeItem::getValue(){
	//return value;
	return text(1).toInt();
}

int IntTreeItem::getMin(){
	return validator->bottom();
}

int IntTreeItem::getMax(){
	return validator->top();
}*/

void IntTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name", text(0));
	xmlWriter.writeAttribute("value", text(1));
	xmlWriter.writeAttribute("type", "int");
}

void IntTreeItem::readInput(QXmlStreamReader &xmlReader){
	setText(1, xmlReader.attributes().value("value").toString());
}

