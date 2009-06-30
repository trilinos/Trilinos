#ifndef TIVABUENA_TREEVIEW_HPP_
#define TIVABUENA_TREEVIEW_HPP_
#include <QTreeView>
#include <QQueue>
namespace TivaBuena{

class Delegate;
class TreeModel;

/**
 * Class used to view TreeModels
 */
class TreeView : public QTreeView{
	Q_OBJECT
public:
	/**
	 * A pair representing an invalidIndex and why it's invalid
	 */
	typedef std::pair<QModelIndex, QString> invalidIndex;

	/**
	 * Constructs a TreeView.
	 * 
	 * @param treeModel The Tree Model being used with the TreeView.
	 * @param delegate The delegate to be used with the TreeView.
	 */
	TreeView(TreeModel *treeModel, Delegate *delegate);

public slots:
	/**
	 * Used to change the visiblity of a row from hidden to shown.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be shown.
	 */
	void showRow(int row, const QModelIndex& parent);

	/**
	 * Used to change the visiblity of a row from shown to hidden.
	 *
	 * @param row The row to be shwon.
	 * @param parent The parent of the item to be hidden.
	 */
	void hideRow(int row, const QModelIndex& parent);

	/**
	 * Handles any badValue signals that might be emitted by the
	 * TreeModel.
	 *
	 * @param badValueIndex The index of the item with the bad value.
	 * @param A brief message explaining what happened to cause the
	 * treeitem to have an invalid value.
	 */
	void handleBadValue(QModelIndex badValueIndex, QString message);

	/**
	 * Checks to see if there are any other invalid indicies.
	 * If there are, it dequeues the next invalidIndex from the 
	 * invalidIndicies queue and calls the handleBadValue function
	 * with it.
	 */
	void checkForOtherBadValues();

private:
	/**
	 * A Queue containing any invalid indicies that need to be
	 * delt with.
	 */
	QQueue<invalidIndex> invalidInicies; 
};



}
#endif //TIVABUENA_TREEVIEW_HPP_
