// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef __im_exodus_ext_h
#define __im_exodus_ext_h

#define IM_EX_INQ_EDGE             27              /* inquire number of edges    */
#define IM_EX_INQ_EDGE_BLK         28              /* inquire number of edge     */
                                                /*   blocks                   */
#define IM_EX_INQ_EDGE_SETS        29              /* inquire number of edge     */
                                                /*   sets                     */
#define IM_EX_INQ_ES_LEN           30              /* inquire length of concat   */
                                                /*   edge set edge list       */
#define IM_EX_INQ_ES_DF_LEN        31              /* inquire length of concat   */
                                                /*   edge set dist factor list*/
#define IM_EX_INQ_EDGE_PROP        32              /* inquire number of props    */
                                                /*   stored per edge block    */
#define IM_EX_INQ_ES_PROP          33              /* inquire number of props    */
                                                /*   stored per edge set      */
#define IM_EX_INQ_FACE             34              /* inquire number of faces    */
#define IM_EX_INQ_FACE_BLK         35              /* inquire number of face     */
                                                /*   blocks                   */
#define IM_EX_INQ_FACE_SETS        36              /* inquire number of face     */
                                                /*   sets                     */
#define IM_EX_INQ_FS_LEN           37              /* inquire length of concat   */
                                                /*   face set face list       */
#define IM_EX_INQ_FS_DF_LEN        38              /* inquire length of concat   */
                                                /*   face set dist factor list*/
#define IM_EX_INQ_FACE_PROP        39              /* inquire number of props    */
                                                /*   stored per face block    */
#define IM_EX_INQ_FS_PROP          40              /* inquire number of props    */
                                                /*   stored per face set      */
#define IM_EX_INQ_ELEM_SETS        41              /* inquire number of face     */
                                                /*   sets                     */
#define IM_EX_INQ_ELS_LEN          42              /* inquire length of concat   */
                                                /*   face set face list       */
#define IM_EX_INQ_ELS_DF_LEN       43              /* inquire length of concat   */
                                                /*   face set dist factor list*/
#define IM_EX_INQ_ELS_PROP         44              /* inquire number of props    */
                                                /*   stored per elem set      */
#define IM_EX_INQ_EDGE_MAP         45              /* inquire number of edge     */
                                                /*   maps                     */
#define IM_EX_INQ_FACE_MAP         46              /* inquire number of face     */
                                                /*   maps                     */

  /*   properties               */
#define IM_EX_EDGE_BLOCK           6               /* edge block property code   */
#define IM_EX_EDGE_SET             7               /* edge set property code     */
#define IM_EX_FACE_BLOCK           8               /* face block property code   */
#define IM_EX_FACE_SET             9               /* face set property code     */
#define IM_EX_ELEM_SET            10               /* face set property code     */
#define IM_EX_EDGE_MAP            11               /* edge map property code     */
#define IM_EX_FACE_MAP            12               /* face map property code     */
#define IM_EX_GLOBAL              13               /* global "block" for variables*/
#define IM_EX_NODAL               14               /* nodal "block" for variables*/



#endif /* __exodus_ext_h */
