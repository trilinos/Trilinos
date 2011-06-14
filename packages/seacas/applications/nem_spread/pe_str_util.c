/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#include <stdio.h>
#include <string.h>
#include <ctype.h>

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s):   Gary L. Hennigan (SNL 9221)
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *      token_compare()
 *      strip_string()
 *      clean_string()
 *      string_to_lower()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int token_compare(char *token, const char *key)
{

  int i1, key_len, kcnt=0;

  key_len = strlen(key);

  for(i1=0; i1 < strlen(token); i1++)
  {
    if(isupper(token[i1]))
      token[i1] = tolower(token[i1]);

    if(token[i1] != ' ')
    {
      if(token[i1] == key[kcnt])
      {
        kcnt++;
        if(kcnt > key_len)
          return 0;
      }
      else
        return 0;
    }
    if(key[kcnt] == ' ')
      kcnt++;
  }

  if(kcnt == strlen(key))
    return 1;
  else
    return 0;

} /*--------------End token_compare()-----------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void strip_string(char inp_str[], const char *tokens)
{
  int  i, j, itok, ntokes, bval;

  i = 0;
  ntokes = strlen(tokens);

  while(inp_str[i] != '\0')
  {
    bval = 0;
    for(itok=0; itok < ntokes; itok++)
    {
      if(inp_str[i] == tokens[itok])
      {
        i++;
        bval = 1;
        break; /* out of for loop */
      }
    }
    if(bval == 0)
      break; /* out of while loop */
  }

  /* Move real part of string to the front */
  j = 0;
  while(inp_str[j+i] != '\0')
  {
    inp_str[j] = inp_str[j+i];
    j++;
  }
  inp_str[j] = inp_str[j+i];
  j--;

  /* Remove trailing tokens */
  while(j != -1)
  {
    bval = 0;
    for(itok=0; itok < ntokes; itok++)
    {
      if(inp_str[j] == tokens[itok])
      {
        bval = 1;
        j--;
        break; /* out of for loop */
      }
    }
    if(bval == 0)
      break; /* out of while loop */
  }

  inp_str[j+1] = '\0';

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void clean_string(char inp_str[], const char *tokens)
{
  int  i, j, itok, ntokes, bval, inplen;

  ntokes = strlen(tokens);
  inplen = strlen(inp_str);

  i = 0;
  bval = 0;
  while(inp_str[i] != '\0')
  {
    for(itok=0; itok < ntokes; itok++)
    {
      if(inp_str[i] == tokens[itok])
      {
        /* Find out if the next character is also a token */
        for(j=0; j < ntokes; j++)
        {
          if(inp_str[i+1] == tokens[j])
          {
            bval = 1;
            break;
          }
        }

        if(bval == 1)
        {
          for(j=i+1; j < inplen; j++)
            inp_str[j] = inp_str[j+1];

          inplen--;
          bval = 0;
          i--;
	  if (i < 0) i = 0;
        }
      }
    }

    i++;

  } /* End "while(inp_str[i] != '\0')" */

  return;

} /*---------------- End clean_string() -----------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void string_to_lower(char in_string[], const char cval)
{
  int len, cnt;

  len = strlen(in_string);
  for(cnt=0; cnt < len; cnt++)
  {
    if(in_string[cnt] == cval)
      return;

    if(isupper(in_string[cnt]))
      in_string[cnt] = tolower(in_string[cnt]);
  }

  return;
}

