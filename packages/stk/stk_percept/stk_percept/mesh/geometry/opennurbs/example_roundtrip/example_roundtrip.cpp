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
    {
      model.Polish();
      if ( model.IsValid() )
      {
        dump->Print("Model is valid after calling Polish().\n");
      }
      else
      {
        dump->Print("Model is not valid.\n");
      }
    }

    /*
    int oi = 14;
    if ( oi >=0 && oi < model.m_object_table.Count() )
    {
      dump->Print("m_object_table[%d].m_object:\n",oi);
      dump->PushIndent();
      model.m_object_table[oi].m_object->Dump(*dump);
      dump->PopIndent();
    }
    */

    // create a text dump of the model
    int version = 4; // File can be read by Rhino 4 and Rhino 5
    //int version = 5; // File can be read by Rhino 5

    ON_String outfile = sFileName;
    int len = outfile.Length() - 4;
    outfile.SetLength(len);
    outfile += "_roundtrip.3dm";
    bool outrc = model.Write( outfile, version, "roundtrip", dump );
    if ( outrc )
    {
      dump->Print("model.Write(%s) succeeded.\n",outfile.Array());
      ONX_Model model2;
      if ( model2.Read( outfile, dump ) )
      {
        dump->Print("model2.Read(%s) succeeded.\n",outfile.Array());
        if ( model2.IsValid(dump) )
        {
          dump->Print("Model2 is valid.\n");
        }
        else
        {
          dump->Print("Model2 is not valid.\n");
        }
        /*
        if ( oi >=0 && oi < model2.m_object_table.Count() )
        {
          dump->Print("m_object_table[%d].m_object:\n",oi);
          dump->PushIndent();
          model2.m_object_table[oi].m_object->Dump(*dump);
          dump->PopIndent();
        }
        */
      }
      else
      {
        dump->Print("model2.Read(%s) failed.\n",outfile.Array());
      }
    }
    else
      dump->Print("model.Write(%s) failed.\n",outfile.Array());

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

