/*
 * solveritem.cpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */

#include "solveritem.hpp"
#include "aztecooitem.hpp"
#include "belositem.hpp"
#include "amesositem.hpp"

SolverItem::SolverItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, SolverItem::Type)
{
	setText(0,"Linear Solver Type");
	setText(1,"AztecOO");
	setText(2,"String");

	solverChild = new AztecOOItem(this);

	solversList << "Amesos" << "Belos";
}

QStringList SolverItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Solver";
	return menuOptions;
}

QString SolverItem::getOptionType(QString optionName){
	return "List";
}

QString SolverItem::getOptionLabel(QString optionName){
	return "Solver:";
}

QStringList SolverItem::getOptionOptions(QString optionName){
	return solversList;
}

void SolverItem::changeParameter(QString optionName, QString value){
	if(value != text(1)){
		solversList << text(1);
		if(value == "AztecOO"){
			delete solverChild;
			solverChild = new AztecOOItem(this);
			addChild(solverChild);
		}
		else if(value == "Belos"){
			delete solverChild;
			solverChild = new BelosItem(this);
			addChild(solverChild);
		}
		else if(value == "Amesos"){
			delete solverChild;
			solverChild = new AmesosItem(this);
			addChild(solverChild);
		}
		solversList.removeAll(value);
		setText(1,value);
	}
}
void SolverItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name",text(0));
	xmlWriter.writeAttribute("type",text(2));
	xmlWriter.writeAttribute("value",text(1));
	xmlWriter.writeStartElement("ParameterList");
	xmlWriter.writeAttribute("name", "Linear Solver Types");
	solverChild->writeOutput(xmlWriter);
	xmlWriter.writeEndElement();
}

void SolverItem::readInput(QXmlStreamReader &xmlReader){
	changeParameter("Change Solver", xmlReader.attributes().value("value").toString());	
	ParameterListTreeItem::readInput(xmlReader);
}

