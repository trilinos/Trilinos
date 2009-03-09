/*
 * belossolveritem.cpp
 *
 *  Created on: Dec 11, 2008
 *      Author: kurtis
 */

#include "belossolveritem.hpp"
#include "blockcgitem.hpp"
#include "blockgmresitem.hpp"
#include "psuedoblockgmresitem.hpp"

BelosSolverItem::BelosSolverItem(ParameterTreeItem *parent)
	:ParameterListTreeItem(parent, BelosSolverItem::Type)
{
	setText(0,"Solver Type");
	setText(1,"Block CG");
	setText(2,"string");

	solverChild = new BlockCGItem(this);

	solversList << "Block GMRES" << "Psuedo Block GMRES" << "GCRODR";
}

QStringList BelosSolverItem::getMenuOptions(){
	QStringList menuOptions = QStringList() << "Change Solver";
	return menuOptions;
}

QString BelosSolverItem::getOptionType(QString optionName){
	return "List";
}

QString BelosSolverItem::getOptionLabel(QString optionName){
	return "Solver:";
}

QStringList BelosSolverItem::getOptionOptions(QString optionName){
	return solversList;
}

void BelosSolverItem::changeParameter(QString optionName, QString value){
	if(value != text(1)){
		solversList << text(1);
		if(value == "Block CG"){
			delete solverChild;
			solverChild = new BlockCGItem(this);
			addChild(solverChild);
		}
		else if(value == "Block GMRES"){
			delete solverChild;
			solverChild = new BlockGMRESItem(this);
			addChild(solverChild);
		}
		else if(value == "Psuedo Block GMRES"){
			delete solverChild;
			solverChild = new PsuedoBlockGMRESItem(this);
			addChild(solverChild);
		}
		else if(value == "GCRODR"){
			delete solverChild;
			QList<ParameterTreeItem*> dummyList;
			solverChild = new ParameterListTreeItem("GCRODR",this, dummyList);
			addChild(solverChild);
		}
		solversList.removeAll(value);
		setText(1,value);
	}
}

void BelosSolverItem::writeOutput(QXmlStreamWriter &xmlWriter){
	xmlWriter.writeEmptyElement("Parameter");
	xmlWriter.writeAttribute("name","Solver Type");
	xmlWriter.writeAttribute("type",text(2));
	xmlWriter.writeAttribute("value",text(1));
	xmlWriter.writeStartElement("ParameterList");
	xmlWriter.writeAttribute("name", "Solver Types");
	solverChild->writeOutput(xmlWriter);
	xmlWriter.writeEndElement();
}

void BelosSolverItem::readInput(QXmlStreamReader &xmlReader){
	changeParameter("Change Solver", xmlReader.attributes().value("value").toString());	
	ParameterListTreeItem::readInput(xmlReader);
}

