#include <iostream>
#include <ostream>
#include <sstream>

#include <streambuf>

#include <QTextEdit>
/**
 * Allows the redirection of output from standard out to a QTextEdit
*/
class CustomStreamBuff: public std::basic_streambuf< char, std::char_traits< char > >
{
	typedef std::basic_streambuf< char, std::char_traits< char > > :: int_type int_type;
	typedef std::char_traits< char > traits_t;
public:
	/**
	 * Constructs a CustomStreamBuff object.
	 *
	 * @param outputDisplay The QTextEdit to which the output should go.
	 */
	CustomStreamBuff(QTextEdit *outputDisplay){
		this->outputDisplay = outputDisplay;
	}
private:
	/**
	 * Acutally outputs to the QTextEdit
	 */
	int_type overflow(int_type c){
		std::cerr << traits_t::to_char_type(c);
		outputDisplay->insertPlainText(QString(c));
		return c;
	}

	QTextEdit *outputDisplay;
};
