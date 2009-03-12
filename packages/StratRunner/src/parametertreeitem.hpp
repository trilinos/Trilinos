/*
 * parameter.hpp
 *
 *  Created on: Dec 6, 2008
 *      Author: kurtis
 */

#ifndef PARAMETERTREEITEM_HPP_
#define PARAMETERTREEITEM_HPP_
#include <QTreeWidgetItem>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
class QTreeWidget;
class SolverTree;

/**
 * ParameterTreeItem is a fundamental class in the stratimikosGUI. All items in the SolverTree are of type ParameterTreeItem. ParameterTreeItem is essentially QTreeWidgetItem with add functions to help the SolverTree display context menus. 
 */
class ParameterTreeItem : public QTreeWidgetItem{
public:
	/**
	 * Constructs a ParameterTreeItem object.
	 *
	 * @param parent The SolverTree of which this ParameterTreeItem is a top level item.
	 * @param type The type of ParameterTreeItem
	 */
	ParameterTreeItem(SolverTree *parent, int type);
	/**
	 * Constructs a ParameterTreeItem object.
	 *
	 * @param parent The parent ParameterTreeItem. 
	 * @param type The type of ParameterTreeItem.
	 */
	ParameterTreeItem(ParameterTreeItem *parent, int type);
	/**
	 * Retruns a list of all the options (if any) that should be a displayed in a context menu if this ParameterTreeItem is currently selected.
	 *
	 * @Return  a list of all the options that should be a displayed in a context menu if this ParameterTreeItem is currently selected.
	 */
	virtual QStringList getMenuOptions()=0;
	/**
	 * Returns the type of option specified by the optionName parameter. The optionName should be one of the options specified in the QStringList returned by the getMenuOptions function. Valid return values are:
	 * <ul>
	 * 	<li>int</li>
	 * 	<li>double</li>
	 * 	<li>bool</li>
	 * 	<li>List</li>
	 * 	<li>filename</li>
	 * </ul>
	 *
	 * @param optionName The option for which the type is desired.
	 * @return The type of the option specified by the optionName parameter.
	 */
	virtual QString getOptionType(QString optionName)=0;
	/**
	 * Returns the label (if any) that should be used in the dialog box that would be used for this option. Note, not all options neccessarily will used a dialog box, and thus this fucntion would return an empty string in such a case.
	 *
	 * @param optionName The name of the option for which the label is desired.
	 * @return The label to be used in the dialog box for the option specified by the optionName parameter.
	 */
	virtual QString getOptionLabel(QString optionName)=0;
	/**
	 * This function retreives any extraneous information about an option that might be important to the option being displayed in a dialog box. For instance, if the option is of type "List" then this function will return a list of all the valid choices to be used for that parameter.
	 *
	 * @param optionName The name of the option for which information is whished to be known.
	 * @return Extreteous but important information about the option specified in the optionName parameter that is important for its display in a dialog box.
	 */
	virtual QStringList getOptionOptions(QString optionName)=0;
	/**
	 * Changes the value of the specified option.
	 *
	 * @param option The option whose value should be changed.
	 * @param The value to which the option should now be set.
	 */
	virtual void changeParameter(QString option, QString value)=0;
	/**
	 * Writes the output for this ParameterTreeItem to the output file.
	 *
	 * @param writer The output stream to which the ParameterTreeItem should write.
	 */
	virtual void writeOutput(QXmlStreamWriter &xmlWriter)=0;
	/**
	 * Reads in a saved .xml file and sets the parameter (or parameters) for the ParameterTreeItem to their correct values based on the input file.
	 *
	 * @param xmlReader The input stream from which to read the file.
	 */
	virtual void readInput(QXmlStreamReader &xmlReader)=0;
protected:
	/**
	 * Finds the child item of this ParameterTreeItem with the specified name. Return 0 if no child is found with the name specified in the arguments.
	 *
	 * @param name The name of the child for which you are looking.
	 * @return The child ParameterTreeItem with the name matching the argument in the name parameter.
	 */
	ParameterTreeItem* findChildByName(QString name);
};

#endif /* PARAMETERTREEITEM_HPP_ */
