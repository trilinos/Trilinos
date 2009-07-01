#ifndef METAWINDOW_HPP_
#define METAWINDOW_HPP_
#include <QMainWindow>
#include <QDialog>
#include <QModelIndex>
#include "Teuchos_ParameterList.hpp"

class QAction;
class QMenu;
class QLabel;
class QPushButton;
class QLineEdit;
namespace TivaBuena{

class TreeModel;
class Delegate;
class TreeView;
class DependencySheet;
/**
 * A small widget that searchs through a parameter list for a particular name
 * of either a parameter or another parameter list.
 */
class SearchWidget : public QDialog{
	Q_OBJECT
public:
	/**
	 * Constructs a SearchWidget.
	 *
	 * @param treeModel The TreeModel being searched.
	 * @param treeView The TreeView being used to display the model.
	 * @param parent The parent widget.
	 */
	SearchWidget(TreeModel *treeModel, TreeView *treeView, QWidget *parent=0);

private slots:
	/**
	 * Searches the for a parameter or parameter list containing the string enterd
	 * in the search terms box.
	 */
	void search();

	/**
	 * Highlights the next result in the list of results that are set
	 * by the search function.
	 */
	void next();

	/**
	 * Highlights the previous result in the list of results that are set
	 * by the search function.
	 */
	void previous();

private:
	/**
	 * Widgets comprising a search widget
	 */
	QPushButton *searchButton, *closeButton, *nextButton, *previousButton;
	QLineEdit *searchTermsEdit;
	QLabel *matchesLabel;
	TreeModel *treeModel;
	TreeView *treeView;

	/**
	 * The results of the search last performed.
	 */
	QList<QModelIndex> currentSearchResults;

	/**
	 * An iterator over the results of the last search performed.
	 */
	QList<QModelIndex>::const_iterator currentSearchIterator;
};

/**
 * The Main Window that contains all other widgets in the stratimikosGUI.
 * For all undocumented functions please refer to the Qt API.
 */
class MetaWindow : public QMainWindow{
	Q_OBJECT
public:
	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param fileName The name of a save file that may store previous values used by a user for the Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	QString fileName=QString());

	/**
	 * Constructs a MainWindow object.
	 * 
	 * @param validParameters The Parameter List the metawindow will display and the user will edit.
	 * @param dependencySheet A sheet listing any dependencies between parameters in the validParameters
	 * ParameterList.
	 * @param fileName The name of a save file that may store previous values used by a user for the Parameter List specified by validParameters.
	 */
	MetaWindow(Teuchos::RCP<Teuchos::ParameterList> validParameters, 
	Teuchos::RCP<DependencySheet> dependencySheet,
	QString fileName=QString());

	/**
	 * Deconstructer for the metawindow
	 */
	~MetaWindow();

protected:
	/**
	 * Handles any QCloseEvents for the metawindow.
	 *
	 * @param event The QCloseEvent that was issued.
	 */
	void closeEvent(QCloseEvent *event);


private:
	/**
	 * Widgets comprising the MetaWindow
	 */
	SearchWidget *searchWidget;
	QAction *resetAct, *loadAct, *saveAct, *saveAsAct, *quitAct, *aboutAct, *searchAct;
	QMenu *fileMenu, *recentMenu, *helpMenu;

	/**
	 * Load and save directory paths
	 */
	QString currentLoadDir, currentSaveDir;

	/**
	 * A list of recently used documents.
	 */
	QStringList recentDocsList;

	/**
	 * The TreeView being used in the metawindow.
	 */
	TreeView *view;

	/**
	 * The TreeModel being used to display the inputs.
	 */
	TreeModel *model;

	/**
	 * The deleages being used to modify any input values.
	 */
	Delegate *delegate;

	/**
	 * Common initialization shared by both constructors
	 */
	void initilization();

	/**
	 * Creates all the menus for the metawindow.
	 */
	void createMenus();

	/**
	 * Creates all necessary actions used in the menut items.
	 */
	void createActions();

	/**
	 * Loads previous parameter settings
	 */
	void load();

	/**
	 * Loads the last state of the MetaWindow (things like window size and screen position).
	 */
	void loadLastSettings();

	/**
	 * Saves the state of the MetaWindow (things like window size and screen position).
	 */
	void saveSettings();

	/**
	 * Currently under developement
	 */
	void addRecentDocument(QString recentDocument);

	/**
	 * Currently under developement
	 */
	void updateRecentDocsMenu();

private slots:
	/**
	 * Creates a new stratimikos solver
	 */
	void resetModel();

	/**
	 * Saves the current solver settings to a user specified file.
	 */
	bool saveFileAs();

	/**
	 * Saves the current solver to the file the user has already specified.
	 */
	void saveFile();

	/**
	 * Loads a solver the user was previously working on and had saved.
	 */
	void loadFile();

	/**
	 * Asks the user whether or not they would like to currently save the file they are working on.
	 * Should be used when the user has modified the file and is about to perform an action that would cause those modifiation to be lost.
	 */
	bool saveCurrentUnsavedFile();

	/**
	 * Loads a document from the set of recent documents
	 */
	void loadRecentDoc();

	/**
	 * Shows information about the program.
	 */
	void showAbout();

	/**
	 * Starts a search for a parituclar Parameter or ParameterList.
	 */
	void initiateSearch();
};



}
#endif /* METAWINDOW_HPP_ */
