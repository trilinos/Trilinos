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
    
    $ python rol_install.py ${{INSTALL_PATH}} ${{OPTS}}
    
    If INSTALL_PATH is not specified, the script will 
    attempt to install to
    
    /usr/local/include 
   
    Install options can be set through OPTS:

    Build Option         OPTS argument            Default Implementation
    --------------------------------------------------------------------
    ROL::Ptr             'shared_ptr'             Teuchos::RCP
    ROL::ParameterList   'property_tree'          Teuchos::ParameterList
    ROL::stacktrace      'backward_cpp'           Teuchos::stacktrace
                         'no_stacktrace'
    
    """

    print(notice)
 
    # Number of command line arguments
    narg = len(sys.argv)
 
    install_specified = ( narg >= 2 )

    src_path = os.path.join(os.getcwd(),"src")
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
        opt_text += "\t Enabled options: \n"

    # ROL::Ptr Implementation
    shared_ptr    = 'shared_ptr'    in options

    if shared_ptr: 
        opt_text += "\t Use std::shared_ptr as ROL::Ptr\n"
    else:
        opt_text += "\t Use Teuchos::RCP as ROL::Ptr\n"


    # ROL::ParameterList Implementation
    property_tree = 'property_tree' in options

    if property_tree: 
        opt_text += "\t Use boost::property_tree as ROL::ParameterList\n"
    else:
        opt_text += "\t Use Teuchos::ParameterList as ROL::ParameterList\n"

    # ROL::stacktrace Implementation
    backward_cpp  = 'backward_cpp'  in options
    no_stacktrace = 'no_stacktrace' in options

    if backward_cpp:
      opt_text += "\t Use backward-cpp for ROL::stacktrace\n"
    elif no_stacktrace:
      opt_text += "\t ROL::stacktrace disabled\n"
    else:
      opt_text += "\t Use Teuchos::stacktrace for ROL::stacktrace\n"


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
    
    """.format(src_path,install_path)

    print(status) 


    numfiles = 0

    for root, dirs, files in os.walk(src_path):
        path = root.split(os.sep)
        headers = [ os.path.join(root,file) for file in files if '.hpp' in file and file[0] != '.' ]


        # Exclude testproblems
        headers = [ h for h in headers if "testproblems" not in h ]

        if shared_ptr:
            headers = [ h for h in headers if 'rcp' not in h ]
        else:
            headers = [ h for h in headers if 'shared_ptr' not in h ]
    
        if property_tree:
            headers = [ h for h in headers if 'parameterlist' not in h ]
        else:
            headers = [ h for h in headers if 'property_tree' not in h ]

        if backward_cpp:
            headers = [ h for h in headers if 'noop' not in h 
                                          and 'teuchos/stacktrace' not in h ]
        elif no_stacktrace:
            headers = [ h for h in headers if 'backward_cpp' not in h 
                                          and 'teuchos/stacktrace' not in h ]
        else: 
            headers = [ h for h in headers if 'backward_cpp' not in h 
                                          and 'noop' not in h ]
        for h in headers:
            print("Copying {0}".format(os.path.split(h)[-1]))
            copy(h, install_path)
            numfiles += 1
    result = \
    """
    
    Copied {0} ROL header files from: {1}
                                  to: {2}

    {3}
    
    """.format(numfiles,src_path,install_path,opt_text)
    print(result)
    print("\nInstallation successful.\n")

    rol_txt = open("ROL.txt",'r')
    rol_logo = rol_txt.read()
    print(rol_logo)
    
    
