import argparse
import os
import re
from shutil import copyfile, copytree
import subprocess
import sys

try:
    doc_path = f"{os.sep}".join(os.path.abspath(__file__).split(os.sep)[:-2])
    sys.path.append(doc_path)
    project_path = f"{os.sep}".join(os.path.abspath(__file__).split(os.sep)[:-4])
    sys.path.append(project_path)
except Exception as e:
    print(f"Can not add project path to system path! Exiting!\nERROR: {e}")
    exit(1)


def is_rst_file(file_path: str) -> bool:
    """ Checks if file_path has .rst extension and if file_path is a file.
    """
    if os.path.splitext(file_path)[-1] == '.rst' and os.path.isfile(file_path):
        return True
    return False


def change_paths_and_get_includes(source_file: str, src_file_path: str, start_path: str, rst_dir: str,
                                  tribits_base_dir: str, copy_file: bool = True) -> tuple:
    """ Changes paths in source file, to be relative to sphinx_path or parent .rst document.
        Returns a tuple with .rst file content and includes(absolute_path, relative_to).
    """
    with open(source_file, 'r') as src_file:
        source_file_str = src_file.read()
        source_file_list = list()
        include_file_list = set()
        for line in source_file_str.split('\n'):
            splitted_line = line.split()
            if 'include::' in splitted_line:
                incl_index = splitted_line.index('include::')
                path_index = incl_index + 1
                if len(splitted_line) > path_index:
                    new_line = []
                    spaces_indented = line.find('.')  # Below 'join()' statement adds a space!
                    if spaces_indented > 0:
                        new_line.append(' ' * (spaces_indented - 1))
                    new_line.extend(splitted_line[:path_index])
                    abs_path = os.path.abspath(os.path.join(src_file_path,
                                                            splitted_line[path_index]))
                    tbd = tribits_base_dir.split(os.sep)[1:]
                    path_elem = abs_path.split(os.sep)[len(tbd) + 1:]
                    new_path = os.path.join(rst_dir, *path_elem)
                    os.makedirs(os.path.dirname(new_path), exist_ok=True)
                    if not os.path.isfile(new_path) and copy_file:
                        copyfile(src=abs_path, dst=new_path, follow_symlinks=True)
                    if is_rst_file(file_path=new_path):
                        include_file_list.add(abs_path)
                    rel_path_from_sphinx_dir = os.path.relpath(path=new_path, start=start_path)
                    new_line.append(rel_path_from_sphinx_dir)
                    new_line = ' '.join(new_line)
                    # Make sure version is formatted correctly:
                    if ':Version:' in new_line:
                        new_line = f":Version:\n  {new_line.split(':Version:')[-1]}"
                    source_file_list.append(new_line)
                else:
                    source_file_list.append(line)
            else:
                source_file_list.append(line)
        abs_path_str = '\n'.join(source_file_list)

    return abs_path_str, include_file_list


class SphinxRstGenerator:
    """ Changes include paths to relative to Sphinx build dir. Saves three main .rst docs files inside Sphinx dir.
    """

    def __init__(self):
        self.paths = {
            'maintainers_guide': {
                'src': os.path.join(doc_path, 'guides', 'maintainers_guide', 'TribitsMaintainersGuide.rst'),
                'src_path': os.path.join(doc_path, 'guides', 'maintainers_guide'),
                'final_path': os.path.join(doc_path, 'sphinx', 'maintainers_guide', 'index.rst'),
                'sphinx_path': os.path.join(doc_path, 'sphinx', 'maintainers_guide'),
                'title': 'TriBITS Maintainers Guide and Reference'},
            'users_guide': {
                'src': os.path.join(doc_path, 'guides', 'users_guide', 'TribitsUsersGuide.rst'),
                'src_path': os.path.join(doc_path, 'guides', 'users_guide'),
                'final_path': os.path.join(doc_path, 'sphinx', 'users_guide', 'index.rst'),
                'sphinx_path': os.path.join(doc_path, 'sphinx', 'users_guide'),
                'title': 'TriBITS Users Guide and Reference'},
            'build_ref': {
                'src': os.path.join(doc_path, 'build_ref', 'TribitsBuildReference.rst'),
                'src_path': os.path.join(doc_path, 'build_ref'),
                'final_path': os.path.join(doc_path, 'sphinx', 'build_ref', 'index.rst'),
                'sphinx_path': os.path.join(doc_path, 'sphinx', 'build_ref'),
                'title': 'Generic TriBITS Project, Build, Test, and Install Reference Guide'}}
        self.rst_dir = os.path.join(doc_path, 'sphinx', 'copied_files')
        self.tribits_base_dir = self._cli()
        self.already_modified_files = set()
        self.create_rst_dir()
        self.build_docs()

    @staticmethod
    def _cli() -> str:
        """ Support for common line arguments. """
        parser = argparse.ArgumentParser()
        parser.add_argument("--copy-base-dir", help="Path to TriBITS base directory")
        args = parser.parse_args()
        abs_path = os.path.abspath(args.copy_base_dir)
        if not abs_path or not os.path.exists(abs_path):
            print(f"\n==> Path: `{abs_path}` is not correct!")
            sys.exit(1)
        print(f"Provided TriBITS base dir: {abs_path}")
        return abs_path

    def create_rst_dir(self) -> None:
        """ Creates copied_files directory in Sphinx directory. All include files will be copy
          there.
        """
        if self.rst_dir is not None:
            if not os.path.exists(self.rst_dir):
                os.makedirs(self.rst_dir)

    @staticmethod
    def build_docs() -> None:
        """ Builds TriBITS documentation based on shell scripts.
        """
        build_script_path = os.path.join(doc_path, 'build_docs.sh')
        subprocess.call([build_script_path, '--skip-final-generation'])

    @staticmethod
    def run_sphinx(cwd: str) -> None:
        """ Runs Sphinx for each documentation.
        """
        sphinx_command = ["make", "html"]
        subprocess.call(sphinx_command, cwd=cwd)

    def combine_documentation(self, docs_dir: str, change_url_to_landing_page: bool = True,
                              change_title_of_docs_main_page: bool = True, title: str = '') -> None:
        """ Renames and moves directory of generated static pages into combined directory
        """
        new_name = os.path.split(docs_dir)[-1]
        dir_to_rename = os.path.join(docs_dir, '_build', 'html')
        new_name_path = os.path.join(docs_dir, '_build', new_name)
        os.rename(src=dir_to_rename, dst=new_name_path)
        static_dir = os.path.join(doc_path, 'sphinx', 'combined_docs', new_name)
        copytree(src=new_name_path, dst=static_dir)
        if change_url_to_landing_page:
            self.change_url_to_landing_page(docs_static_dir=static_dir)
        if change_title_of_docs_main_page:
            self.change_title_of_docs_main_page(docs_static_dir=static_dir, new_title=title)

    @staticmethod
    def change_url_to_landing_page(docs_static_dir: str) -> None:
        """ Changes home url of documentation page, so it points to landing page.
        """
        index_html = os.path.join(docs_static_dir, 'index.html')
        with open(index_html, 'r') as index_read:
            index_str = index_read.read()
        repl_url = index_str.replace('<a href="#" class="icon icon-home"> TriBITS',
                                     '<a href="../index.html" class="icon icon-home"> TriBITS')
        with open(index_html, 'w') as index_write:
            index_write.write(repl_url)

    @staticmethod
    def change_title_of_docs_main_page(docs_static_dir: str, new_title: str) -> None:
        """ Changes home url of documentation page, so it points to landing page.
        """
        index_html = os.path.join(docs_static_dir, 'index.html')
        with open(index_html, 'r') as index_read:
            index_str = index_read.read()
        repl_url = re.sub('(<title>1 Introduction &mdash;).+(documentation</title>)',
                          f'<title>{new_title}</title>', index_str)
        with open(index_html, 'w') as index_write:
            index_write.write(repl_url)

    @staticmethod
    def save_rst(file_path: str, file_content: str) -> None:
        """ Saves .rst file with given pathh and content
        """
        with open(file_path, 'w') as dest_file:
            dest_file.write(file_content)

    def generate_rst(self, source_file: str, final_path: str = None, src_path: str = None,
                     start_path: str = None) -> set:
        """ Generate correct links in .rst files, so Sphinx can find them
        """
        if final_path is None:
            overwrite_source = True
        else:
            overwrite_source = False

        file_content, includes = change_paths_and_get_includes(source_file=source_file, src_file_path=src_path,
                                                               start_path=start_path, rst_dir=self.rst_dir,
                                                               tribits_base_dir=self.tribits_base_dir)

        if overwrite_source:
            self.save_rst(file_path=source_file, file_content=file_content)
        else:
            self.save_rst(file_path=final_path, file_content=file_content)

        return includes

    def remove_title_numbering(self) -> None:
        """ Removes numbering from docs.
        """
        for doc_name, sources in self.paths.items():

            str_to_replace = '.. rubric::'
            with open(sources.get('final_path'), 'r') as src_file:
                org_str = src_file.read()
                org_list = org_str.split('\n')
                if org_list[0].startswith('====='):
                    del org_list[0]
                if org_list[1].startswith('====='):
                    del org_list[1]
                org_list[0] = f'{str_to_replace} {org_list[0]}'
                mod_str = '\n'.join(org_list)

            with open(sources.get('final_path'), 'w') as dst_file:
                dst_file.write(mod_str)

    def main(self):
        """ Main routine goes for nested .rst docs
        """
        child_rst = set()
        for doc_name, sources in self.paths.items():
            includes = self.generate_rst(source_file=sources.get('src'), src_path=sources.get('src_path'),
                                         final_path=sources.get('final_path'), start_path=sources.get('sphinx_path'))
            child_rst.update(includes)
        self.already_modified_files.update(child_rst)
        tbd = self.tribits_base_dir.split(os.sep)[1:]
        child_rst_lst = list(child_rst)

        sphinx_rel_path = self.paths.get('maintainers_guide').get('sphinx_path')
        grand_child_rst = set()
        for child in child_rst_lst:
            path_elem = child.split(os.sep)[len(tbd) + 1:]
            final_path = os.path.join(self.rst_dir, *path_elem)
            os.makedirs(os.path.dirname(final_path), exist_ok=True)
            src_path = os.path.split(child)[0]
            includes_grand = self.generate_rst(source_file=child, src_path=src_path,
                                               final_path=final_path, start_path=sphinx_rel_path)
            grand_child_rst.update(includes_grand)
        grand_child_rst_lst = [gc_rst for gc_rst in grand_child_rst if gc_rst not in self.already_modified_files]

        grand_grand_child_rst = set()
        for grand_child in grand_child_rst_lst:
            path_elem = grand_child.split(os.sep)[len(tbd) + 1:]
            final_path = os.path.join(self.rst_dir, *path_elem)
            os.makedirs(os.path.dirname(final_path), exist_ok=True)
            src_path = os.path.split(grand_child)[0]
            includes_grand_grand = self.generate_rst(source_file=grand_child, src_path=src_path,
                                                     final_path=final_path, start_path=sphinx_rel_path)
            grand_grand_child_rst.update(includes_grand_grand)

        if not grand_grand_child_rst:
            print('DONE! ALL GOOD!\n')
        else:
            print('NOT DONE!\n')

        self.remove_title_numbering()

        print('===> Generating Sphinx documentation:\n')
        for doc_name, sources in self.paths.items():
            cwd = sources.get('sphinx_path')
            print(f'===> Generating {doc_name}\n')
            self.run_sphinx(cwd=cwd)
            self.combine_documentation(docs_dir=cwd, title=sources.get('title'))


if __name__ == '__main__':
    SphinxRstGenerator().main()
