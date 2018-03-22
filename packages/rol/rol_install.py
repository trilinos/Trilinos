#!/usr/bin/python

import sys, os, errno

from shutil import copy

if __name__ == '__main__':

    notice = \
    """
    ***************************************************
    *** Installing Rapid Optimization Library (ROL) ***
    ***        Header-only installation             ***
    ***************************************************
    
    Example Usage:
    
    $ python rol_install.py ${{INSTALL_PATH}} OPTS
    
    If INSTALL_PATH is not specified, the script will 
    attempt to install to
    
    /usr/local/include 
    
    Optional arguments can also be provided through OPTS.
    Currently the only supported option is shared_ptr 
    to use the std::shared_ptr instead of Teuchos::RCP.
    
    """

    print(notice)
 
    # Number of command line arguments
    narg = len(sys.argv)
 
    install_specified = ( narg >= 2 )

    src_path = os.getcwd()
    install_path = '/usr/local/include'
    options = sys.argv[2:]

    if narg > 1:
        install_path = sys.argv[1]

    # Check if install path exists
    if not os.path.exists( install_path ): 
        print("\nInstall path {0} does not exist. Attempting to create it.".format(install_path))
   
        try:
            os.makedirs(install_path)
        except OSError as e:
            if e.errno == errno.EACCES:
                print("\nYou do not have necessary permission to create path {0}!".format(install_path))
                print("\nROL installation failed! Exiting.\n")
                os._exit(1)
            elif e.errno != errno.EEXIST:
                print("\nUnable to create path {0}. Try another path.".format(install_path))
                print("\nROL installation failed! Exiting.\n")
                os._exit(1)
            pass
        print("\nPath {0} created.".format(install_path))    
    
    install_fail_msg = \
    """
    Install of ROL failed! 
    You do not have write permissions to the designated install path: 
    
    {0}\n""".format(install_path)

    if not os.access(install_path, os.W_OK):
        print( install_fail_msg ) 
        os._exit(1)

    opt_text = ""

    if len(options)>0:
        opt_text += "Enabled options: \n"

    shared_ptr = 'shared_ptr' in options

    if shared_ptr: 
        opt_text += "-> Use std::shared_ptr as ROL::Ptr"

    status = \
    """
    Main source directory (where all ROL directories live):
    ----> {0}
    
    Source /src directory (where we'll get headers from):
    ----> {0}/src
    
    Install directory (the main installation directory):
    ----> {1}
    
    Include directory (where we'll install headers):
    ----> {1}/include
    
    {2}
    """.format(src_path,install_path,opt_text)

    print(status) 


    numfiles = 0

    for root, dirs, files in os.walk(src_path):
        path = root.split(os.sep)
        headers = [ os.path.join(root,file) for file in files if '.hpp' in file and file[0] != '.' ]

        if shared_ptr:
            headers = [ h for h in headers if 'rcp' not in h ]
        else:
            headers = [ h for h in headers if 'shared_ptr' not in h ]
    
        for h in headers:
            copy(h, install_path)
            print("Copying {0}".format(os.path.split(h)[-1]))
            numfiles += 1
    result = \
    """
    
    Copied {0} ROL header files from: {1}
                                  to: {2}
    
    """.format(numfiles,src_path,install_path)
    print("\nInstallation successful.\n")

    rol_txt = open("ROL.txt",'r')
    rol_logo = rol_txt.read()
    print(rol_logo)
    
    
