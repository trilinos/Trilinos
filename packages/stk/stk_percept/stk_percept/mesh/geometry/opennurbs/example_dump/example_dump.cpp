// uncomment the "ON_DLL_IMPORTS" define to use opennurbs as a Windows DLL
//#define ON_DLL_IMPORTS
#include "../opennurbs.h"
#include "../examples_linking_pragmas.h"

ON_BOOL32 Dump3dmFile( 
        const char*, // full name of file
        ON_TextLog&
        );

static ON_BOOL32 bTerseReport = 0;

int main ( int argc, const char* argv[] )
{
  ON::Begin();

  // default dump is to stdout
  // use the -out:filename.txt option to dump to a file
  ON_TextLog dump_to_stdout;
  ON_TextLog* dump = &dump_to_stdout;
  FILE* dump_fp = 0;
  
  dump->SetIndentSize(2);

  int argi;
  if ( argc < 2 ) 
  {
    printf("Syntax: %s [-out:outputfilename.txt] [-terse] file1.3dm file2.3dm ...\n",argv[0]);
    return 0;
  }

  for ( argi = 1; argi < argc; argi++ ) 
  {
    const char* arg = argv[argi];
    if ( 0 == arg )
      continue;

    // check for -terse option
    if ( 0 == strcmp( arg, "-terse" ) ) 
    {
      bTerseReport = 1;
      continue;
    }

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
          dump = 0;
          ON::CloseFile(dump_fp);
          dump_fp = 0;
        }
        dump_fp = text_fp;
        text_fp = 0;
        dump = new ON_TextLog(dump_fp);
      }
      continue;
    }

    Dump3dmFile( arg, *dump );
    dump->Print("\n\n");
  }

  if ( dump_fp )
  {
    // close the text dump file
    delete dump;
    dump = 0;
    ON::CloseFile( dump_fp );
    dump_fp = 0;
  }

  // OPTIONAL: Call just before your application exits to clean
  //           up opennurbs class definition information.
  //           Opennurbs will not work correctly after ON::End()
  //           is called.
  ON::End();

  return 0;
}

ON_BOOL32 Dump3dmFile( const char* sFileName, ON_TextLog& dump )
{
  dump.Print("====== FILENAME: %s\n",sFileName);
  ON_Workspace ws;
  FILE* fp = ws.OpenFile( sFileName, "rb" ); // file automatically closed by ~ON_Workspace()
  if ( !fp ) {
    dump.Print("**ERROR** Unable to open file.\n");
    return false;
  }

  ON_BinaryFile file( ON::read3dm, fp );

  int version = 0;
  ON_String comment_block;
  ON_BOOL32 rc = file.Read3dmStartSection( &version, comment_block );
  if (!rc) {
    dump.Print("**ERROR** Read3dmStartSection() failed\n");
    return false;
  }
  dump.Print("====== VERSION: %d\n",version );
  dump.Print("====== COMMENT BLOCK:\n",version );
  dump.PushIndent();
  dump.Print(comment_block);
  dump.PopIndent();
  dump.Print("====== CHUNKS:\n",version );
  unsigned int typecode;
  while ( !file.AtEnd() ) {
    typecode = file.Dump3dmChunk( dump, 0 );
    if ( !typecode )
      break;
    if ( typecode == TCODE_ENDOFFILE )
      break;
  }
  dump.Print("====== FINISHED: %s\n",sFileName);

  return true;
}

