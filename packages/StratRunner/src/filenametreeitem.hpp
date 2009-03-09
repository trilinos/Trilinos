/*
 * stringedittreeitem.hpp
 *
 *  Created on: Dec 19, 2008
 *      Author: kurtis
 */

#ifndef FILENAMETREEITME_HPP_
#define FILENAMETREEITME_HPP_
#include "parametertreeitem.hpp"

/**
 * FileNameTreeItem is a ParameterTreeItem that represents a file name parameter used in Stratimikos
 */
class FileNameTreeItem : public ParameterTreeItem{
public:
	enum {Type = UserType+19};
	/**
	 * Constructs a FileNameTreeItem object.
	 *
	 * @param name The name of the file name parameter (will be displayed to the user).
	 * @param dialogLabel The label that will be displayed next to the input region in the dialog that pops up when
	 * the user attemps to change the value of the parameter.
	 * @param parent The parent ParameterTreeItem.
	 */
	FileNameTreeItem(QString name, QString dialogLabel, ParameterTreeItem *parent);
	QStringList getMenuOptions();
	QString getOptionType(QString optionName);
	QString getOptionLabel(QString optionName);
	QStringList getOptionOptions(QString optionName);
	void changeParameter(QString option, QString value);
	void changeValue(QString value);
	void writeOutput(QXmlStreamWriter &xmlWriter);
	void readInput(QXmlStreamReader &xmlReader);
private:
	QString name, dialogLabel;
};

#endif /* FILENAMETREEITME_HPP_ */

