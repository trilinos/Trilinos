#ifndef TIVABUENA_TREEITEM_HPP_
#define TIVABUENA_TREEITEM_HPP_

#include <QList>
#include <QVariant>
#include <QXmlStreamWriter>
#include <QXmlStreamReader>
#include "Teuchos_ParameterList.hpp"
#include "TivaBuena_ArrayHelperFunctions.hpp"
namespace TivaBuena{
/**
 * The TreeItem class is the item class used by the TreeModel class.
 */
class TreeItem{
public:
	/**
	 * Constructs a TreeItem object.
	 *
	 * @param data A list of data that should be in the TreeItem. The list should be of length 3 and contain the following data in 
	 * each respective location:
	 * <ol>
	 *	<li>The name of the parameter</li>
	 *	<li>The default value of the parameter</li>
	 *	<li>The type of parameter</li>
	 * </ol>
	 * In the case of a TreeItem representing a ParameterList the data list should contain the following in each
	 * respective location:
	 * <ol>
	 * 	<li>The name of the ParameterList</li>
	 * 	<li>An empty string</li>
	 * 	<li>The "list" parameter type</li>
	 * </ol>
	 * @param parameterEntry The ParameterEntry this TreeItem is ment to represent.
	 * @param parent The parent TreeItem.
	 * @param unrecognized If true, this item will be unrecognized and not displayed, if false the item will be displayed.
	 */
	TreeItem(const QList<QVariant> &data, Teuchos::ParameterEntry *parameterEntry, TreeItem *parent = 0, bool unrecognized=false);

	/**
	 * Deconstrcutor for the TreeItem.
	 */
	~TreeItem();

	/**
	 * Prints out the values in the TreeItem.
	 */
	void printOut() const;

	/**
	 * Appends a child TreeItem to the TreeItem
	 * 
	 * @param child The child item to be appended.
	 */
	void appendChild(TreeItem *child);

	/**
	 * Returns the child treeitem in the row specified by the row argument.
	 *
	 * @param row The row in which the child is in.
	 * @return The child TreeItem located in the row.
	 */
	TreeItem *child(int row);

	/**
	 * Gets the number of child nodes this item has.
	 *
	 * @return The number of child nodes this item has.
	 */
	int childCount() const;

	/**
	 * Gets a list of all the child items.
	 *
	 * @return A list of all child items.
	 */
	const QList<TreeItem*> getChildItems();

	/**
	 * How man columns the TreeItem has. Should always be 3.
	 *
	 * @return The number of columns the TreeItem has.
	 */
	int columnCount() const;

	/**
	 * Returns the data located in a particular column.
	 *
	 * @param column The column of the desired data.
	 * @param role The role of the data.
	 * @return The data located in a particular column.
	 */
	QVariant data(int column, int role = Qt::DisplayRole) const;

	/**
	 * Gets the parent TreeItem
	 *
	 * @return The parent TreeItem.
	 */
	TreeItem *parent();

	/**
	 * Returns the row in which this TreeItem is located.
	 * 
	 * @return The row in which this TreeItem is located.
	 */
	int row() const;

	/**
	 * Gets the ParameterEntry associated with this TreeItem if it has one.
	 *
	 * @return The ParameterEntry associated with this TreeItem if it has one.
	 */
	const Teuchos::ParameterEntry *entry();

	/**
	 * Determines whether or not the current value associated with the
	 * TreeItem is valid.
	 *
	 * @return True if the value is valid, false otherwise.
	 */
	bool hasValidValue() const;

	/**
	 * Writes the current state of the TreeItem to an xml stream
	 *
	 * @param xmlWriter The xml stream to which the object should write.
	 */
	void writeOutput(QXmlStreamWriter &xmlWriter);

	/**
	 * Changes the value of the TreeItem. Should only be used with TreeItems that represent Parameters.
	 *
	 * @param value The new value to be assigned to the TreeItem.
	 */
	bool changeValue(QVariant value);

	/**
	 * Sets the validator for the parameter the TreeItem represents.
	 *
	 * @param validator The validator which the parameter should be given.
	 */
	void setValidator(Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator);

private:
	/**
	 * Changes the value of an array.
	 *
	 * @param value A string representing the value of the array.
	 * @param type The type of array.
	 */
	void changeValueForArray(QVariant value, QString type);

	/**
	 * Whether or not the parameter type is unrecognized.
	 */
	bool unrecognized;

	/**
	 * The childitems of the TreeItem.
	 */
	QList<TreeItem*> childItems;

	/**
	 * The data in the item.
	 */
	QList<QVariant> itemData;

	/**
	 * The parent TreeItem.
	 */
	TreeItem *parentItem;

	/**
	 * The ParameterEntry being represented by the TreeItem.
	 */
	Teuchos::ParameterEntry *parameterEntry;

	/**
	 * The docString for the TreeItem.
	 */
	QString docString;
};



}
#endif /* TIVABUENA_TREEITEM_HPP_ */
