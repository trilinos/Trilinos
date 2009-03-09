/*
 * parameerlisttreeitem.hpp
 *
 *  Created on: Dec 10, 2008
 *      Author: kurtis
 */

#ifndef PARAMEERLISTTREEITEM_HPP_
#define PARAMEERLISTTREEITEM_HPP_
#include "parametertreeitem.hpp"
#include <QList>
class SolverTree;

/**
 * ParameterListTreeItem is a fundamental class in the stratimikosGUI. ParameterListTreeItems have a "blank" context menu. They simply contain other ParameterTreeItems and other ParameterListTreeItems.
 */
class ParameterListTreeItem : public ParameterTreeItem{
public:
	/**
	 * Constructs a ParameterTreeItem object.
	 *
	 * @param parent The parent ParameterTreeItem. 
	 * @param type The type of ParameterTreeItem.
	 */
	ParameterListTreeItem(ParameterTreeItem *parent, int type);
	/**
	 * Constructs a ParameterListTreeItem object.
	 *
	 * @param parent The SolverTree of which this ParameterListTreeItem is a top level item.
	 * @param type The type of ParameterListTreeItem
	 */
	ParameterListTreeItem(SolverTree *parent, int type);
	/**
	 * Constructs a ParameterTreeItem object.
	 *
	 * @param parent The parent ParameterTreeItem. 
	 * @param type The type of ParameterTreeItem.
	 * @param children A list of ParameterTreeItems that should be initially added to the ParameterListTreeItem.
	 */
	ParameterListTreeItem(QString name, ParameterTreeItem *parent, QList<ParameterTreeItem*> children);
	/**
	 * Constructs a ParameterListTreeItem object.
	 *
	 * @param parent The SolverTree of which this ParameterListTreeItem is a top level item.
	 * @param type The type of ParameterListTreeItem
	 * @param children A list of ParameterTreeItems that should be initially added to the ParameterListTreeItem.
	 */
	ParameterListTreeItem(QString name, SolverTree *parent, QList<ParameterTreeItem*> children);
	/**
	 * Returns a blank list of options since a ParameterListTreeItem does not have changeable options.
	 *
	 * @return A blank list of Menu Options
	 */
	virtual QStringList getMenuOptions();
	/**
	 * Returns a blank string since there are not options for a ParameterListTreeItem
	 *
	 * @return A blank sting.
	 */
	virtual QString getOptionType(QString optionName);
	/**
	 * Returns a blank string since there are not options for a ParameterListTreeItem
	 *
	 * @return A blank string.
	 */
	virtual QString getOptionLabel(QString optionName);
	/**
	 * Returns a blank list since a ParameterListTreeItem does not have changeable options.
	 *
	 * @return A blank list.
	 */
	virtual QStringList getOptionOptions(QString optionName);
	/**
	 * Does nothing since there are no options to change.
	 */
	virtual void changeParameter(QString optionName, QString value);
	virtual void writeOutput(QXmlStreamWriter &xmlWriter);
	virtual void readInput(QXmlStreamReader &xmlReader);
protected:
	/**
	 * Adds all the parameters that should be children of this ParameterListTreeItem. By default, no parameters are add to the ParameterListTreeItem. 
	 *
	 * If ParameterListTreeItem is being subclassed, this function should be redefined in the sublcass and then called in the constructor of the subclass.
	 */
	virtual void addParameters();
	/**
	 * Adds all the parameter lists that should be children of this ParameterListTreeItem. By default, no parameter lists are add to the ParameterListTreeItem. 
	 *
	 * If ParameterListTreeItem is being subclassed, this function should be redefined in the sublcass and then called in the constructor of the subclass.
	 */
	virtual void addParameterLists();
private:
	/**
	 * Adds items to the ParameterListTreeItem that are specified by the children parameter.
	 *
	 * @param children Children to be added to the ParameterListTreeItem.
	 */
	void addChildren(QList<ParameterTreeItem*> children);
};

#endif /* PARAMEERLISTTREEITEM_HPP_ */
