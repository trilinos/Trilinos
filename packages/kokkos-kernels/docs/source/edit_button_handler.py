import re
import os
import sys

try:
    project_path = f'{os.sep}'.join(os.path.abspath(__file__).split(os.sep)[:-1])
    print(f'==> Starting in: {project_path}')
    sys.path.append(project_path)
except Exception as e:
    print(f'Can not add project path to system path! Exiting!\nERROR: {e}')
    raise SystemExit(1)


class FileFinder:
    """Finds files with given extension in given directory and subdirectories as well. Returns a list of them."""
    def __init__(self, directory: str, file_extension: list):
        self.directory = os.path.abspath(directory)
        self.file_extension = file_extension
        self.files_list = None

    def get_files(self) -> list:
        """Returning found files."""
        if self.files_list is None:
            self.files_list = []
            for path, subdirs, files in os.walk(self.directory):
                for name in files:
                    if os.path.isfile(os.path.join(path, name)) and os.path.splitext(os.path.join(path, name))[-1] \
                            in self.file_extension:
                        self.files_list.append(os.path.join(path, name))
        return self.files_list


class HTMLButtonAdder:
    """Adds button for direct edition docs pages on GitHub."""
    def __init__(self, document_files: list, html_files: list, excluded_files: list, html_tag: str):
        self.document_files = sorted(document_files)
        self.html_files = sorted(html_files)
        self.excluded_files = excluded_files
        self.html_tag = html_tag
        self.__initial_checks()

    def __initial_checks(self):
        """Performs initial checks. Makes sure adding an edit button run smoothly."""
        print('==> Starting initial checks:')
        # Making sure the document_files has the same length as html_files
        if len(self.html_files) != len(self.document_files):
            raise AssertionError(f'Length of `html_files` list: {len(self.html_files)} is different than length of '
                                 f'`document_files` list: {len(self.document_files)}. Must be the same. '
                                 f'Check excluded files.')
        print(f'=> Found {len(self.html_files)} files')
        # Making sure that string to replace could be found in each html_file
        html_tag_list = []
        missing_files = []
        for file in self.html_files:
            with open(file) as html_file:
                html_str = html_file.read()
                if html_str.find(self.html_tag) != -1:
                    html_tag_list.append(file)
                else:
                    missing_files.append(file)
        if len(self.html_files) != len(html_tag_list):
            raise AssertionError(f'String to replace was not found in files: {missing_files}')
        print(f'=> Found {len(html_tag_list)} files with matching string to replace')
        print('--------------------------------------------------')

    def __overwrite_html(self, file_names: tuple, wiki_prefix: str) -> None:
        """Overwriting html file with button addition."""
        # Reading file, replacing string and overwriting
        with open(file_names[0], 'rt') as html_file:
            html_file_str = html_file.read()
        # Preparing icon string
        replaced_str = file_names[1].replace(project_path, wiki_prefix)
        str_to_put = f'\n<div class="edit-this-page">\n  <a class="muted-link" href="{replaced_str}" title="Edit ' \
                     f'this page">\n    <svg aria-hidden="true" viewBox="0 0 24 24" stroke-width="1.5" ' \
                     f'stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">\n      ' \
                     f'<path stroke="none" d="M0 0h24v24H0z" fill="none"></path>\n      <path d="M4 20h4l10.5' \
                     f' -10.5a1.5 1.5 0 0 0 -4 -4l-10.5 10.5v4"></path>\n      <line x1="13.5" y1="6.5" x2="17.5" ' \
                     f'y2="10.5"></line>\n    </svg>\n    <span class="visually-hidden">Edit this page</span>\n  ' \
                     f'</a>\n</div>\n'
        if html_file_str.find(str_to_put) != -1:
            print(f'=> File: {file_names[0]} already has an edit button')
            return
        # Finding position to put edit icon
        tag_pos = html_file_str.find(self.html_tag)
        html_str_beg = html_file_str[:tag_pos + len(html_tag)]
        html_str_end = html_file_str[tag_pos + len(html_tag):]
        html_str_replace = html_str_beg + str_to_put + html_str_end
        html_str_replace = self._overwrite_deprecated_style(html_string=html_str_replace)
        with open(file_names[0], 'wt') as new_html_file:
            new_html_file.write(html_str_replace)
        print(f'=> Processing: {file_names[0]} done')

    @staticmethod
    def _overwrite_deprecated_style(html_string: str) -> str:
        """Overwrites deprecated style."""
        # Overwriting the deprecated style without version
        str_to_replace = '<span class="pre">[DEPRECATED]</span>'
        replaced_with = '<span class="pre" style="color:#A020F0;font-weight:bold;">[DEPRECATED]</span>'
        html_string = html_string.replace(str_to_replace, replaced_with)
        # Overwriting the deprecated style with version
        not_finished = True
        deprecated_regex = '<span class="pre">\[DEPRECATED<\/span> <span class="pre">since<\/span> ' \
                           '<span class="pre">([0-9]+.[0-9]+.*[0-9]*){1}]<\/span>'
        while not_finished:
            match = re.search(deprecated_regex, html_string, re.MULTILINE)
            if match is None:
                return html_string
            else:
                html_str_beg = html_string[:match.regs[0][0]]
                html_str_end = html_string[match.regs[0][1]:]
                html_tag = html_string[match.regs[0][0]:match.regs[0][1]]
                html_tag = html_tag.replace('"pre"', '"pre" style="color:#A020F0;font-weight:bold;"')
                html_string = html_str_beg + html_tag + html_str_end


    def add_button(self, wiki_prefix: str) -> None:
        """Loops over html files and overwrite them."""
        for num, file_names in enumerate(zip(self.html_files, self.document_files), 1):
            print(f'==> Processing pair {num}:\n=> {file_names[0]}\n=> {file_names[1]}')
            self.__overwrite_html(file_names=file_names, wiki_prefix=wiki_prefix)


if __name__ == "__main__":
    print('==================================================')
    print('==> Starting adding buttons to html files:')
    # Getting lists of documents and html files
    document_files = FileFinder(directory=project_path, file_extension=['.md', '.rst']).get_files()
    generated_docs_dir = os.path.join(project_path, '../generated_docs/docs')
    html_files = FileFinder(directory=generated_docs_dir, file_extension=['.html']).get_files()
    # Excluded files (Files created by Sphinx, not to be overwritten with edit button)
    excluded_files = [os.path.abspath(os.path.join(generated_docs_dir, 'genindex.html')),
                      os.path.abspath(os.path.join(generated_docs_dir, 'search.html'))]
    # Final `html_files` list of files to add edit button to
    html_files = [html_file for html_file in html_files if html_file not in excluded_files]
    # HTML tag after edit button is injected
    html_tag = '<div class="content-icon-container">'
    print(f'=> Adding button after: {html_tag}')
    # Wiki prefix pointing directly to GitHub
    wiki_prefix = 'https://github.com/kokkos/kokkos-core-wiki/blob/main/docs/source'
    print(f'=> Using prefix for Kokkos Wiki: {wiki_prefix}')
    HTMLButtonAdder(document_files=document_files, html_files=html_files, excluded_files=excluded_files,
                    html_tag=html_tag).add_button(wiki_prefix=wiki_prefix)
    print('==================================================')
