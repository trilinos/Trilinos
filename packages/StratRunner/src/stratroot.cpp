/*
 * stratroot.cpp
 *
 *  Created on: Dec 5, 2008
 *      Author: kurtis
 */
#include "stratroot.hpp"
#include "solveritem.hpp"
#include "preconitem.hpp"
#include "solvertree.hpp"
#include "booltreeitem.hpp"
#include <QInputDialog>

StratRoot::StratRoot(SolverTree *parent)
	:ParameterListTreeItem(parent, StratRoot::Type)
{
	setText(0,QString("Stratimikos Solution"));
	setText(1,QString(""));
	setText(2,QString(""));

	addParameters();
	addParameterLists();
	setExpanded(true);
}

void StratRoot::addParameters(){
	addChild(new BoolTreeItem("Enable Delayed Solver Construction", this));
}

void StratRoot::addParameterLists(){
	addChild(new SolverItem(this));
	addChild(new PreconItem(this));
}

void StratRoot::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeStartElement("ParameterList");
	for(int i=0; i<childCount(); i++){
		dynamic_cast<ParameterTreeItem*>(child(i))->writeOutput(xmlWriter);
	}
	xmlWriter.writeEndElement();
}

void StratRoot::readInput(QXmlStreamReader &xmlReader){
	std::string token;
	while(!xmlReader.isEndDocument()){
		token = xmlReader.tokenString().toStdString();	
		if(xmlReader.isStartElement()){
			std::string name = xmlReader.attributes().value("name").toString().toStdString();
			ParameterTreeItem *foundItem = 	findChildByName(xmlReader.attributes().value("name").toString());
			if(foundItem != 0){
				foundItem->readInput(xmlReader);
			}
		}
		xmlReader.readNext();
	}
}

