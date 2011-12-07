/* $Header: /src4/opennurbs/example_read/example_read.cpp 2     8/29/06 1:43p Dalelear $ */
/* $NoKeywords: $ */
//
// Copyright (c) 1993-1998 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
////////////////////////////////////////////////////////////////
//
//  example_read.cpp  
// 
//  Example program using the Rhino file IO toolkit.  The program reads in  
//  a Rhino 3dm model file and describes its contents.  The program is a 
//  console application that takes a filename as a command line argument.
//
////////////////////////////////////////////////////////////////////////

// uncomment the "ON_DLL_IMPORTS" define to use opennurbs as a Windows DLL
//#define ON_DLL_IMPORTS
#include "../opennurbs.h"
#include "../examples_linking_pragmas.h"

int main( int argc, const char *argv[] )
{
  // If you are using OpenNURBS as a Windows DLL, then you MUST use
  // ON::OpenFile() to open the file.  If you are not using OpenNURBS
  // as a Windows DLL, then you may use either ON::OpenFile() or fopen()
  // to open the file.

  int argi;
  if ( argc < 2 ) 
  {
    printf("Syntax: %s [-out:outputfilename.txt] file1.3dm file2.3dm ...\n",argv[0] );
    return 0;
  }

  // Call once in your application to initialze opennurbs library
  ON::Begin();

  // default dump is to stdout
  ON_TextLog dump_to_stdout;
  ON_TextLog* dump = &dump_to_stdout;
  FILE* dump_fp = 0;

  ONX_Model model;

  bool bVerboseTextDump = true;

  for ( argi = 1; argi < argc; argi++ ) 
  {
    const char* arg = argv[argi];

    // check for -out or /out option
    if ( ( 0 == strncmp(arg,"-out:",5) || 0 == strncmp(arg,"/out:",5) ) 
         && arg[5] )
    {
      // change destination of dump file
      const char* sDumpFilename = arg+5;
      FILE* text_fp = ON::OpenFile(sDumpFilename,"w");
      if ( text_fp )
      {
        if ( dump_fp )
        {
          delete dump;
          ON::CloseFile(dump_fp);
        }
        dump_fp = text_fp;
        dump = new ON_TextLog(dump_fp);
      }
      continue;
    }

    const char* sFileName = arg;

    dump->Print("\nOpenNURBS Archive File:  %s\n", sFileName );

    // open file containing opennurbs archive
    FILE* archive_fp = ON::OpenFile( sFileName, "rb");
    if ( !archive_fp ) 
    {
      dump->Print("  Unable to open file.\n" );
      continue;
    }

    dump->PushIndent();

    // create achive object from file pointer
    ON_BinaryFile archive( ON::read3dm, archive_fp );

    // read the contents of the file into "model"
    bool rc = model.Read( archive, dump );

    // close the file
    ON::CloseFile( archive_fp );

    // print diagnostic
    if ( rc )
      dump->Print("Successfully read.\n");
    else
      dump->Print("Errors during reading.\n");

    // see if everything is in good shape
    if ( model.IsValid(dump) )
      dump->Print("Model is valid.\n");
    else
      dump->Print("Model is not valid.\n");

    // create a text dump of the model
    if ( bVerboseTextDump )
    {
      dump->PushIndent();
      model.Dump(*dump);
      dump->PopIndent();
    }

    // destroy this model
    model.Destroy();

    dump->PopIndent();
  }

  if ( dump_fp )
  {
    // close the text dump file
    delete dump;
    ON::CloseFile( dump_fp );
  }
  
  // OPTIONAL: Call just before your application exits to clean
  //           up opennurbs class definition information.
  //           Opennurbs will not work correctly after ON::End()
  //           is called.
  ON::End();

  return 0;
}

