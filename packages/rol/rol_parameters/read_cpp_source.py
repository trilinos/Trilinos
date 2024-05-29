

import pathlib


def contains_escaped_quote_advanced( s : str ) -> bool:
    """
    Determines if a string contains an escaped double quote character.

    This function checks for occurrences of double quotes (") that are
    preceded by an odd number of backslashes (\), indicating that the
    quote is escaped.

    Parameters:
    s (str): The input string to check.

    Returns:
    bool: True if an escaped double quote is found, False otherwise.
    """
    i = 0
    while i < len(s):
        if s[i] == '\\':
            backslash_count = 1
            i += 1

            # Count consecutive backslashes
            while i < len(s) and s[i] == '\\':
                backslash_count += 1
                i += 1

            # If there's an odd number of backslashes followed by a quote, then it is escaped
            if i < len(s) and s[i] == '"' and backslash_count % 2 == 1:
                return True
        else:
            i += 1
    return False



def strip_cpp_comments( cpp_source : str ) -> str:
    """
    Removes C++ style comments (both single-line and multi-line) from a string of C++ source code.

    This function strips out both single-line (//) and multi-line (/* ... */) comments
    from the provided C++ source code, while preserving the content within string literals.

    Parameters:
    cpp_source (str): The input C++ source code as a string.

    Returns:
    str: The source code with comments removed.
    """    
    in_string = False
    in_single_line_comment = False
    in_multi_line_comment = False
    result = []
    i = 0
    while i < len(cpp_source):
        # Check for string start/end
        if cpp_source[i] == '"' and not (in_single_line_comment or in_multi_line_comment):
            # Extract substring from the current position backwards to the last non-escaped quote or start
            substring = cpp_source[:i+1][::-1]
            # Check if the quote is escaped
            if not contains_escaped_quote_advanced(substring):
                in_string = not in_string
            result.append(cpp_source[i])
        # Check for single-line comment start
        elif i+1 < len(cpp_source) and cpp_source[i:i+2] == "//" and not (in_string or in_multi_line_comment):
            in_single_line_comment = True
            i += 1  # Skip next character to avoid parsing '/' twice
        # Check for multi-line comment start
        elif i + 1 < len(cpp_source) and cpp_source[i:i+2] == "/*" and not (in_string or in_single_line_comment):
            in_multi_line_comment = True
            i += 1  # Skip next character to avoid parsing '*' twice
        # Check for single-line comment end
        elif in_single_line_comment and cpp_source[i] == "\n":
            in_single_line_comment = False
            result.append(cpp_source[i])  # Include newline in result
        # Check for multi-line comment end
        elif i + 1 < len(cpp_source) and in_multi_line_comment and cpp_source[i:i+2] == "*/":
            in_multi_line_comment = False
            i += 1  # Skip next character to avoid parsing '/' twice
        # Append character if not in a comment
        elif not (in_single_line_comment or in_multi_line_comment):
            result.append(cpp_source[i])
        i += 1

    return ''.join(result)



def read_cpp_source( cpp_file : pathlib.Path ) -> str:
    """
    Reads a C++ source file, removes comments, and returns the cleaned source code.

    This function reads the content of a given C++ source file, strips out all comments,
    and returns the resulting cleaned source code as a string.

    Parameters:
    cpp_file (pathlib.Path): The path to the C++ source file to read.

    Returns:
    str: The C++ source code with comments removed.

    Raises:
    AssertionError: If the provided path does not exist or is not a file.
    """
    # Ensure the argument is a file
    assert( cpp_file.exists() )
    assert( cpp_file.is_file() )

    # Read C++ source file to string
    with open(cpp_file,"r") as f:
         content = f.read()

    cpp_source = strip_cpp_comments(content)

    return cpp_source


