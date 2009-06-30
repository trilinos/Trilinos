#ifndef TIVABUEAN_FILENAMEWIDGET_HPP_
#define TIVABUEAN_FILENAMEWIDGET_HPP_
#include <QWidget>
class QLabel;

namespace TivaBuena{

/**
 * A small widget consisting of a button and label that allows the user
 *  to select a file through a QFileDialog. The label displays
 *  the currently selected file.
 */
class FilenameWidget : public QWidget{
	Q_OBJECT
public:
	/**
	 * Constructs a FilenameWidget
	 *
	 * @param currentFileName The Filename with which the widget should be 
	 * initialized.
	 * @param parent The parent widget.
	 */
	FilenameWidget(QString currentFileName=QString(), QWidget* parent=0);

	/**
	 * Gets the current filename in the widget.
	 *
	 * @return The current filename in the widget.
	 */
	QString getCurrentFileName();

	/**
	 * Sets the current filename in the widget.
	 *
	 * @param The name to which the widget should be set.
	 */
	void setCurrentFileName(QString newName);

public slots:
	/**
	 * Opens a QFileDialog allowing the user to select a new filename.
	 */
	void getNewFileName();

private:
	/**
	 * The current file name stored in the list.
	 */
	QString currentFileName;
	/**
	 * The label describing the file path.
	 */
	QLabel *pathLabel;

};


}

#endif //TIVABUEAN_FILENAMEWIDGET_HPP_
