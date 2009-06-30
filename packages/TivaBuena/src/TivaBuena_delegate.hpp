#ifndef TIVABUENA_DELEGATE_HPP_
#define TIVABUENA_DELEGATE_HPP_
#include <QItemDelegate>
#include <QModelIndex>
#include <QSize>
#include "Teuchos_ParameterEntry.hpp"

namespace TivaBuena{

/**
 * The delegate used for the TivaBuena package. For non-documented functions please refer to the Qt API.
 */
class Delegate : public QItemDelegate{
	Q_OBJECT
public:
	/**
	 * Constructs a Delegate.
	 * 
	 * @param parent The parent object of the Delegate.
	 */
	Delegate(QObject *parent = 0);

	QWidget* createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const;
	virtual void setEditorData(QWidget *editor, const QModelIndex &index) const;
	void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const;
	virtual void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const;

private:
	/**
	 * Handles any array editing that needs to be done.
	 *
	 * @param index The index of the item being edited.
	 * @param type The type of array being edited.
	 * @param parent The parent widget.
	 */
	void arrayHandler(const QModelIndex& index, QString type, QWidget *parent) const;
};


}
#endif /* TIVABUENA_DELEGATE_HPP_ */
