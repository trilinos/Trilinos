#include <QtGui>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QComboBox>
#include <QFileDialog>
#include "TivaBuena_delegate.hpp"
#include "TivaBuena_ArrayWidget.hpp"
#include "float.h"
#include <iostream>

namespace TivaBuena{

Delegate::Delegate(QObject *parent):QItemDelegate(parent){}

QWidget* Delegate::createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/ , const QModelIndex &index ) const{
	QWidget *editor = 0;
	if(index.column() != 1)
		return editor;

	Teuchos::RCP<const Teuchos::ParameterEntryValidator> paramValidator = ((TreeModel*)(index.model()))->getValidator(index);
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);

	if(itemType == intId){
		editor = new QSpinBox(parent);
		Teuchos::RCP<const EnhancedNumberValidator <int> > intValidator;
		if(!Teuchos::is_null(paramValidator))
			intValidator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<int> >(paramValidator);
		EnhancedNumberValidator<int>::applyToSpinBox(intValidator, (QSpinBox*)editor);
	}
	else if(itemType == shortId){
		editor = new QSpinBox(parent);
		Teuchos::RCP<const EnhancedNumberValidator<short> > shortValidator;
		if(!Teuchos::is_null(paramValidator))
			shortValidator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<short> >(paramValidator);
		EnhancedNumberValidator<short>::applyToSpinBox(shortValidator, (QSpinBox*)editor);
	}
/*	else if(itemType == longlongId){
		editor = new QwwLongSpinBox(parent);
		Teuchos::RCP<const EnhancedNumberValidator<long long> > longlongValidator;
		if(!Teuchos::is_null(paramValidator))
			longlongValidator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<long long> >(paramValidator);
		EnhancedNumberValidator<long long>::applyToSpinBox(longlongValidator, (QDoubleSpinBox*)editor);
	}*/
	else if(itemType == doubleId){
		editor = new QDoubleSpinBox(parent);
		Teuchos::RCP<const EnhancedNumberValidator<double> > doubleValidator;
		if(!Teuchos::is_null(paramValidator))
			doubleValidator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<double> >(paramValidator);
		EnhancedNumberValidator<double>::applyToSpinBox(doubleValidator, (QDoubleSpinBox*)editor);
	}
	else if(itemType == floatId){
		editor = new QDoubleSpinBox(parent);
		Teuchos::RCP<const EnhancedNumberValidator<float> > floatValidator; 
		if(!Teuchos::is_null(paramValidator))
			floatValidator = Teuchos::rcp_dynamic_cast<const EnhancedNumberValidator<float> >(paramValidator);
		EnhancedNumberValidator<float>::applyToSpinBox(floatValidator, (QDoubleSpinBox*)editor);
	}
	else if(itemType == boolId){
		editor = new QComboBox(parent);
		static_cast<QComboBox*>(editor)->addItem("true");
		static_cast<QComboBox*>(editor)->addItem("false");
	}
	else if(itemType == stringId){
		if(Teuchos::is_null(paramValidator))
			editor = new QLineEdit(parent);
		else if(!Teuchos::is_null(Teuchos::rcp_dynamic_cast<const FileNameValidator>(paramValidator))){
			QString paramName = 
				((TreeModel*)(index.model()))->data(index.sibling(index.row(), 0),Qt::DisplayRole).toString();
			QString currentPath = ((TreeModel*)(index.model()))->data(index,Qt::DisplayRole).toString();
			if(currentPath.size() == 0)
				currentPath = QDir::homePath();
			QString filename = QFileDialog::getSaveFileName(parent, paramName, currentPath);
			if(filename != "")
				((TreeModel*)(index.model()))->setData(index, filename);
		}
		else if(paramValidator->validStringValues()->size() != 0){
			Teuchos::RCP<const Teuchos::Array<std::string> > options = paramValidator->validStringValues();
			editor = new QComboBox(parent);
			for(Teuchos::Array<std::string>::const_iterator itr = options->begin(); itr != options->end(); ++itr){
				static_cast<QComboBox*>(editor)->addItem(QString::fromStdString(*itr));
			}
		}
		else
			editor = new QLineEdit(parent);
	}
	else if(itemType.contains(arrayId))
		arrayHandler(index, itemType.section(" ", -1), parent);

	return editor;
}

void Delegate::arrayHandler(const QModelIndex& index, QString type, QWidget *parent) const{
	if(type == intId){
		IntArrayWidget array(index, type, parent);
		array.exec();
	}
	else if(type == shortId){
		ShortArrayWidget array(index, type, parent);
		array.exec();
	}
	/*else if(type == longlongId){
		LongLongArrayWidget array(index, type, parent);
		array.exec();
	}*/
	else if(type == doubleId){
		DoubleArrayWidget array(index, type, parent);
		array.exec();
	}
	else if(type == floatId){
		FloatArrayWidget array(index, type, parent);
		array.exec();
	}
	else if(type == stringId){
		StringArrayWidget array(index, type, parent);
		array.exec();
	}
}

void Delegate::setEditorData(QWidget *editor, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	if(itemType == intId){
		int value = index.model()->data(index).toInt();
		static_cast<QSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == shortId){
		short value = index.model()->data(index).toInt();
		static_cast<QSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == doubleId){
		double value = index.model()->data(index).toDouble();
		static_cast<QDoubleSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == floatId){
		float value = index.model()->data(index).toDouble();
		static_cast<QDoubleSpinBox*>(editor)->setValue(value);
	}
	else if(itemType == boolId){
		QString value = index.model()->data(index).toString();
		static_cast<QComboBox*>(editor)->setEditText(value);
	}
	else if(itemType == stringId){
		QString value = index.model()->data(index).toString();
		Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = ((TreeModel*)(index.model()))->getValidator(index);
		if(Teuchos::is_null(validator) || validator->validStringValues()->size()==0)
			static_cast<QLineEdit*>(editor)->setText(value);
	 	else
			static_cast<QComboBox*>(editor)->setEditText(value);
	}
}

void Delegate::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const{
	QString itemType = ((TreeModel*)(index.model()))->itemType(index);
	if(itemType == intId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, spinBox->value(), Qt::EditRole);
	}
	else if(itemType == shortId){
		QSpinBox *spinBox = static_cast<QSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, (short)spinBox->value(), Qt::EditRole);
	}
	else if(itemType == doubleId){
		QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, spinBox->value(), Qt::EditRole);
	}
	else if(itemType == floatId){
		QDoubleSpinBox *spinBox = static_cast<QDoubleSpinBox*>(editor);
		spinBox->interpretText();
		model->setData(index, (float)spinBox->value(), Qt::EditRole);
	}
	else if(itemType == boolId){
		bool value = static_cast<QComboBox*>(editor)->currentText() == "true"; 
		model->setData(index, value, Qt::EditRole);
	}
	else if(itemType == stringId){
		Teuchos::RCP<const Teuchos::ParameterEntryValidator> validator = 
			static_cast<const TreeModel*>(index.model())->getValidator(index);
		QString value;
		if(Teuchos::is_null(validator))
			value = static_cast<QLineEdit*>(editor)->text();
		else
			value = static_cast<QComboBox*>(editor)->currentText(); 
		model->setData(index, value, Qt::EditRole);
	}
}

void Delegate::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &/*index*/) const{
	editor->setGeometry(option.rect);
}


}

