/*
 * parameterlisttreeitem.cpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */
#include "parameterlisttreeitem.hpp"
#include "solvertree.hpp"

ParameterListTreeItem::ParameterListTreeItem(SolverTree *parent, int type):ParameterTreeItem(parent, type){

}
ParameterListTreeItem::ParameterListTreeItem(ParameterTreeItem *parent, int type):ParameterTreeItem(parent, type){

}

ParameterListTreeItem::ParameterListTreeItem(QString name, ParameterTreeItem *parent, QList<ParameterTreeItem*> children)
	:ParameterTreeItem(parent, UserType+100)
{
	addChildren(children);
	setText(0,name);
}

ParameterListTreeItem::ParameterListTreeItem(QString name, SolverTree *parent, QList<ParameterTreeItem*> children)
	:ParameterTreeItem(parent, UserType+100)
{
	addChildren(children);
	setText(0,name);
}

void ParameterListTreeItem::addChildren(QList<ParameterTreeItem*> children){
	QList<ParameterTreeItem*>::const_iterator it = children.constBegin();
	while(it != children.constEnd()){
		(*it)->parent()->removeChild(*it);
		addChild(*it);
		++it;
	}
}

QStringList ParameterListTreeItem::getMenuOptions(){
	QStringList menuOptions;
	return menuOptions;
}

QString ParameterListTreeItem::getOptionType(QString optionName){
	return "";
}

QString ParameterListTreeItem::getOptionLabel(QString optionName){
	return "";
}

QStringList ParameterListTreeItem::getOptionOptions(QString optionName){
	QStringList optionOptions;
	return optionOptions;
}

void ParameterListTreeItem::changeParameter(QString optionName, QString value){}

void ParameterListTreeItem::addParameters(){}

void ParameterListTreeItem::addParameterLists(){}

void ParameterListTreeItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeStartElement("ParameterList");
	xmlWriter.writeAttribute("name", text(0));
	for(int i=0; i<childCount(); i++){
		dynamic_cast<ParameterTreeItem*>(child(i))->writeOutput(xmlWriter);
	}
	xmlWriter.writeEndElement();
}

void ParameterListTreeItem::readInput(QXmlStreamReader &xmlReader){
	while(!(xmlReader.isEndElement() && xmlReader.name().toString() == "ParameterList")){
		if(xmlReader.isStartElement()){
			ParameterTreeItem *foundItem = 	findChildByName(xmlReader.attributes().value("name").toString());
			if(foundItem != 0){
				foundItem->readInput(xmlReader);
			}
		}
		xmlReader.readNext();
	}
}

