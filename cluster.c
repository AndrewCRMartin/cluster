/*************************************************************************

   Program:    
   File:       cluster.c
   
   Version:    V2.1
   Date:       20.12.11
   Function:   Perform cluster analysis
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995-2011
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  20.06.95 Original
   V1.1  26.06.95 Fixed potential off-bottom of array accesses.
   V1.2  29.06.95 Fixed bug in freeing memory
   V1.3  03.07.95 Fixed bug in (original code) setting up iorder array
                  for labelling clusters
   V1.4  04.07.95 Added median printing
   V1.5  06.07.95 Some tidying to median printing and prints word
                  MEDIANS:
   V1.6  28.07.95 Added -l   
   V1.7  02.08.95 Fixed bug in calculation of medians
   V1.8  11.09.95 Modified to use fgetsany() and hence avoid any limit
                  on the length of a record.
                  Gets INF from MAXDOUBLE rather than just 1e20
   V2.0  10.06.96 Input format now allows a label at the end of each 
                  vector. This is appended to the clustering table on
                  output.
   V2.1  20.12.11 Modified for new GetWord()

*************************************************************************/
/* Includes
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <values.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"
#include "bioplib/array.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define INF ((REAL)MAXDOUBLE)
#define UP '|'
#define ACROSS '-'
#define BLANK ' '
#define MAXBUFF  160
#define MAXWORD  32

/*  Map row I and column J of upper half diagonal symmetric matrix 
    onto vector.
*/
#define IOFFSET(n,i,j) (j+(i-1)*n-(i*(i+1))/2)

typedef struct _indata
{
   REAL *data;
   struct _indata *next;
}  INDATA;

typedef struct _label
{
   struct _label *next;
   char *label;
}  LABEL;

/************************************************************************/
/* Globals
*/
extern int errno;

/************************************************************************/
/* Prototypes
*/
BOOL HierClus(int n, int m, int iopt, REAL **data, int *ia, int *ib, 
              REAL *crit);
int **ClusterAssign(FILE *fp, int n, int *ia, int *ib, REAL *crit, 
                    int lev, int *iorder, REAL *critval, int *height,
                    LABEL *LabelList);
BOOL InsertIorder(int *iorder, int lev, int cluster, int parent);
char **ClusterDendogram(FILE *fp,int lev, int *iorder, int *height, 
                        REAL *critval);
REAL **ReadData(FILE *in, int *NVec, int *VecDim, LABEL **pLabelList);
BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                  int Method, BOOL ShowTable, BOOL ShowDendogram,
                  int NClus, int NPrint, LABEL *LabelList);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *Method, BOOL *ShowTable, BOOL *ShowDendogram,
                  int *nclus, int *nprint);
BOOL FindMedian(FILE *fp, int **clusters, REAL **data, int NVec, 
                int VecDim, int NClus, int ClusNum);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for cluster analysis

   21.06.95 Original    By: ACRM
   29.06.95 Frees memory allocated for data array
*/
int main(int argc, char **argv)
{
   FILE  *in  = stdin,
         *out = stdout;
   char  InFile[MAXBUFF],
         OutFile[MAXBUFF];
   int   Method,
         NVec,
         VecDim,
         NClus,
         NPrint = 0;
   BOOL  ShowTable,
         ShowDendogram;
   REAL  **data;
   LABEL *LabelList;

   if(ParseCmdLine(argc, argv, InFile, OutFile, &Method, &ShowTable, 
                   &ShowDendogram, &NClus, &NPrint))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((data = ReadData(in, &NVec, &VecDim, &LabelList))!=NULL)
         {
            if(!ShowClusters(out, data, NVec, VecDim, Method, ShowTable, 
                             ShowDendogram, NClus, NPrint, LabelList))
            {
               fprintf(stderr,"Memory allocation for clustering \
failed\n");
               return(1);
            }
            FreeArray2D((char **)data, NVec, VecDim);
         }
         else
         {
            fprintf(stderr,"Memory allocation for reading data failed\n");
            return(1);
         }
      }
      else
      {
         return(1);
      }
   }
   else
   {
      Usage();
   }
   return(0);
}


/************************************************************************/
/*>REAL **ReadData(FILE *in, int *NVec, int *VecDim, LABEL **pLabelList)
   ---------------------------------------------------------------------
   Reads data into a linked list then creates and fills a 2D array from
   these data, freeing the linked list

   21.06.95 Original    By: ACRM
   11.09.95 Modified to use fgetsany()
   10.06.96 Can read a label appended to the vectors. Stores these labels
            in the LabelList linked list which is returned as the new
            parameter pLabelList
   20.12.11 Modified for new GetWord()
*/
#define LabelList (*pLabelList)
REAL **ReadData(FILE *in, int *NVec, int *VecDim, LABEL **pLabelList)
{
   char   *buffer,
          item[MAXWORD],
          *chp,
          *label;
   int    count;
   INDATA *indata, *p;
   REAL   **data;
   BOOL   GotLabel,
          NoLabels = FALSE;
   LABEL  *l;

   LabelList = NULL;
   *NVec     = 0;
   *VecDim   = 0;
   
   /* Get the first line from the file, ignoring comment lines          */
   do 
   {
      if((buffer = fgetsany(in))==NULL)
         return(NULL);
      
      TERMINATE(buffer);
      KILLLEADSPACES(chp,buffer);
   }
   while(*chp == '\0' || *chp == '!' || *chp == '#');
   
   /* Count the number of data items in the line                        */
   GotLabel = FALSE;
   for(chp=buffer, *VecDim=0; chp!=NULL; (*VecDim)++)
   {
      label = chp;
      chp = GetWord(chp, item, MAXWORD);
      if(item[0] == '#')
      {
         GotLabel = TRUE;
         KILLLEADSPACES(label, label+1);
         break;
      }
   }
   
   /* Initialise linked list to store this line                         */
   INIT(indata, INDATA);
   p=indata;
   if(p==NULL) return(NULL);
   if((p->data = (REAL *)malloc(*VecDim * sizeof(REAL)))==NULL)
   {
      free(p);
      return(NULL);
   }
   
   /* Copy in the data                                                  */
   for(chp=buffer, count=0; chp!=NULL; count++)
   {
      chp = GetWord(chp, item, MAXWORD);
      if(item[0] == '#')
         break;
      p->data[count] = (REAL)atof(item);
   }
   *NVec = 1;

   /* Initialise storage for labels                                     */
   INIT(LabelList, LABEL);
   l=LabelList;
   if(LabelList == NULL)
   {
      fprintf(stderr,"cluster: No memory to store label\n");
      NoLabels = TRUE;
   }
   else
   {
      l->label = NULL;
   }
      
   /* Store the label if there is one                                   */
   if(GotLabel && !NoLabels)
   {
      if((l->label = (char *)malloc((strlen(label)+1)*sizeof(char)))
         ==NULL)
      {
         free(LabelList);
         LabelList = NULL;
         fprintf(stderr,"cluster: No memory to store label\n");
         NoLabels = TRUE;
      }
      else
      {
         strcpy(l->label, label);
      }
   }
   
   free(buffer);

   /* Read and store the rest of the data                               */
   while((buffer=fgetsany(in)) != NULL)
   {
      TERMINATE(buffer);
      KILLLEADSPACES(chp, buffer);
      if(*chp == '\0' || *chp == '!' || *chp == '#')
         continue;
      
      (*NVec)++;
      
      /* Allocate and check next position in linked list                */
      ALLOCNEXT(p, INDATA);
      if(p==NULL)
      {
         /* Free the data list                                          */
         for(p=indata; p!=NULL; NEXT(p))
            free(p->data);
         FREELIST(indata, INDATA);
 
         /* Free the label list as well                                 */
         for(l=LabelList; l!=NULL; NEXT(l))
         {
            if(l->label != NULL)
               free(l->label);
         }
         FREELIST(LabelList, LABEL);
         
         return(NULL);
      }

      /* Allocate and check next position in the label list             */
      if(!NoLabels)
      {
         ALLOCNEXT(l, LABEL);
         if(l==NULL)
         {
            /* Free the labels and print a warning                      */
            for(l=LabelList; l!=NULL; NEXT(l))
            {
               if(l->label != NULL)
                  free(l->label);
            }
            FREELIST(LabelList, LABEL);
            LabelList = NULL;

            fprintf(stderr,"cluster: No memory to store label\n");
            NoLabels = TRUE;
         }
         else
         {
            l->label = NULL;
         }
      }

      /* Allocate and check data for this item in linked list           */
      if((p->data = (REAL *)malloc(*VecDim * sizeof(REAL)))==NULL)
      {
         for(p=indata; p!=NULL; NEXT(p))
            free(p->data);
         FREELIST(indata, INDATA);

         /* Free the label list as well                                 */
         for(l=LabelList; l!=NULL; NEXT(l))
         {
            if(l->label != NULL)
               free(l->label);
         }
         FREELIST(LabelList, LABEL);
         LabelList = NULL;
         
         return(NULL);
      }
      
      /* Copy in the data                                               */
      GotLabel = FALSE;
      for(chp=buffer, count=0; chp!=NULL && count<(*VecDim); count++)
      {
         label = chp;
         chp = GetWord(chp, item, MAXWORD);
         if(item[0] == '#') /* This should never happen!                */
         {
            GotLabel = TRUE;
            KILLLEADSPACES(label, label+1);
            break;
         }
         p->data[count] = (REAL)atof(item);
      }

      if(!GotLabel && (chp != NULL))
      {
         label = chp;
         chp = GetWord(chp, item, MAXWORD);
         if(item[0] == '#')
         {
            GotLabel = TRUE;
            KILLLEADSPACES(label, label+1);
         }
      }

      /* Store the label if there is one                                */
      if(GotLabel && !NoLabels)
      {
         if((l->label = (char *)malloc((strlen(label)+1)*sizeof(char)))
            ==NULL)
         {
            /* Free the current label list                              */
            for(l=LabelList; l!=NULL; NEXT(l))
            {
               if(l->label != NULL)
                  free(l->label);
            }
            FREELIST(LabelList, LABEL);
            LabelList = NULL;

            fprintf(stderr,"cluster: No memory to store label\n");
            NoLabels = TRUE;
         }
         else
         {
            strcpy(l->label, label);
         }
      }

      /* Free memory allocated for string                               */
      free(buffer);
   }

   if(errno)  /* Returned by fgetsany()                                 */
   {
      for(p=indata; p!=NULL; NEXT(p))
      {
         free(p->data);
      }
      FREELIST(indata, INDATA);
      perror("System error in cluster");
      return(NULL);
   }
   
   /*** Convert the linked list into a 2D array                       ***/

   /* Allocate and check 2D array                                       */
   if((data = (REAL **)Array2D(sizeof(REAL),*NVec,*VecDim))==NULL)
   {
      for(p=indata; p!=NULL; NEXT(p))
      {
         free(p->data);
      }
      FREELIST(indata, INDATA);
      return(NULL);
   }

   /* Walk the linked list, copying in the data                         */
   for(p=indata,*NVec=0; p!=NULL; NEXT(p),(*NVec)++)
   {
      for(count=0; count<(*VecDim); count++)
         data[*NVec][count] = p->data[count];
      free(p->data);
   }

   /* Free up the linked list (we've already done the actual data)      */
   FREELIST(indata, INDATA);

   /* If there were no labels stored in the label linked list, we'll
      free that as well. If NoLabels is TRUE, then LabelList is
      already NULL.
   */
   if(!NoLabels)
   {
      NoLabels = TRUE;
      for(l=LabelList; l!=NULL; NEXT(l))
      {
         if(l->label != NULL)
         {
            NoLabels = FALSE;
            break;
         }
      }

      /* If no labels were found, then free the list and set to NULL    */
      if(NoLabels)
      {
         FREELIST(LabelList, LABEL);
         LabelList = NULL;
      }
   }
   
   return(data);
}
#undef LabelList


/************************************************************************/
/*>BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                     int Method, BOOL ShowTable, BOOL ShowDendogram,
                     int NClus, int NPrint, LABEL *LabelList)
   -----------------------------------------------------------------
   Allocate temporary arrays, call clustering code and free the arrays

   21.06.95 Original    By: ACRM
   04.07.95 Added NClus and printing of medians
   06.07.95 Prints the word MEDIANS and checks return from FindMedians()
   28.07.95 Added NPrint
   10.06.96 Added LabelList (passed to ClusterAssign())
*/
BOOL ShowClusters(FILE *fp, REAL **data, int NVec, int VecDim, 
                  int Method, BOOL ShowTable, BOOL ShowDendogram, 
                  int NClus, int NPrint, LABEL *LabelList)
{
   int  i,
        *ia,
        *ib,
        lev = (NPrint==0)?NVec:NPrint,
        *iorder,
        *height,
        **clusters;
   REAL *crit,
        *critval;
   BOOL ok = TRUE;
   char **out = NULL;

   ia      = (int *)malloc(NVec * sizeof(int));
   ib      = (int *)malloc(NVec * sizeof(int));
   crit    = (REAL *)malloc(NVec * sizeof(REAL));

   iorder  = (int *)malloc(lev * sizeof(int));
   height  = (int *)malloc(lev * sizeof(int));
   critval = (REAL *)malloc(lev * sizeof(REAL));

   if(ia==NULL      ||
      ib==NULL      ||
      crit==NULL    ||
      critval==NULL ||
      iorder==NULL  ||
      height==NULL)
   {
      if(ia      != NULL) free(ia);
      if(ib      != NULL) free(ib);
      if(iorder  != NULL) free(iorder);
      if(height  != NULL) free(height);
      if(crit    != NULL) free(crit);
      if(critval != NULL) free(critval);
      return(FALSE);
   }
   
   /* Do the clustering                                                 */
   if(HierClus(NVec,VecDim,Method,data,ia,ib,crit))
   {
      /* Assign data to clusters                                        */
      if((clusters = ClusterAssign((ShowTable?fp:NULL),NVec,ia,ib,crit,
                                   lev,iorder,critval,
                                   height,LabelList))!=NULL)
      {
         /* Print the dendogram                                         */
         if(ShowDendogram)
         {
            if((out = ClusterDendogram(fp,lev,iorder,height,critval))
               ==NULL)
               ok = FALSE;
         }

         /* Print the medians                                           */
         if(NClus)
         {
            fprintf(fp,"\nMEDIANS:\n");
            for(i=1; i<=NClus; i++)
            {
               if(!FindMedian(fp, clusters, data, NVec, VecDim, NClus, i))
               {
                  ok = FALSE;
               }
            }
         }
         
      }
      else
      {
         ok = FALSE;
      }
   }
   else
   {
      ok = FALSE;
   }

   /* Free allocated memory                                             */
   if(clusters != NULL) FreeArray2D((char **)clusters,NVec,VecDim);
   if(out      != NULL) FreeArray2D((char **)out,lev*3,lev*3);
   if(ia       != NULL) free(ia);
   if(ib       != NULL) free(ib);
   if(iorder   != NULL) free(iorder);
   if(height   != NULL) free(height);
   if(crit     != NULL) free(crit);
   if(critval  != NULL) free(critval);
   
   return(ok);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     int *Method, BOOL *ShowTable, BOOL *ShowDendogram,
                     int *nclus, int *nlimit)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            REAL   *radius      Neighbour radius
            int    *CountType   Counting scheme
            int    *nclus       Number of clusters for which to show
                                medians
            int    *nprint      Number of levels to display table and
                                dendogram
   Returns: BOOL                Success?

   Parse the command line
   
   20.06.95 Original    By: ACRM
   28.07.95 Added -l and nprint
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *Method, BOOL *ShowTable, BOOL *ShowDendogram,
                  int *nclus, int *nprint)
{
   argc--;
   argv++;

   *Method        = 1;
   *ShowTable     = TRUE;
   *ShowDendogram = TRUE;
   *nclus         = 0;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'm':
            if(--argc == 0)
               return(FALSE);
            argv++;

            if(!upstrncmp(argv[0],"WAR",3) || argv[0][0] == '1')
               *Method = 1;
            else if(!upstrncmp(argv[0],"SIN",3) || argv[0][0] == '2')
               *Method = 2;
            else if(!upstrncmp(argv[0],"COM",3) || argv[0][0] == '3')
               *Method = 3;
            else if(!upstrncmp(argv[0],"AVE",3) || 
                    !upstrncmp(argv[0],"GRO",3) || argv[0][0] == '4')
               *Method = 4;
            else if(!upstrncmp(argv[0],"MCQ",3) || argv[0][0] == '5')
               *Method = 5;
            else if(!upstrncmp(argv[0],"MED",3) || 
                    !upstrncmp(argv[0],"GOW",3) || argv[0][0] == '6')
               *Method = 6;
            else if(!upstrncmp(argv[0],"CEN",3) || argv[0][0] == '7')
               *Method = 7;
            else
               return(FALSE);
            break;
         case 't':
            *ShowTable = FALSE;
            break;
         case 'd':
            *ShowDendogram = FALSE;
            break;
         case 'n':
            argc--;
            argv++;
            if((argc < 0) || (sscanf(argv[0],"%d",nclus) == 0))
               return(FALSE);
            break;
         case 'l':
            argc--;
            argv++;
            if((argc < 0) || (sscanf(argv[0],"%d",nprint) == 0))
               return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL HierClus(int NVec, int VecDim, int ClusterMethod, REAL **data, 
                 int *ia, int *ib, REAL *crit)
   -------------------------------------------------------------------
   Input:   int  NVec                 Number of vectors to cluster
            int  VecDim               Dimension of each vector
            REAL data[NVec][VecDim]   Input data matrix
            int  ClusterMethod        Clustering criterion to be used
   Output:  int  ia[NVec]             \
                 ib[NVec]             | History of allomerations
            REAL crit[NVec]           /
   Returns: BOOL                      Success of memory allocations

   Hierarchical clustering using user-specified criterion. 
                                                             
   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.
   26.06.95 Fixed frees on error
   29.06.95 Fix the pointers in the data array when finished
*/
BOOL HierClus(int NVec, int VecDim, int ClusterMethod, REAL **data, 
              int *ia, int *ib, REAL *crit)
{
   int  ind, 
        ind1, 
        ind2, 
        ind3, 
        NClusters, 
        i, 
        j, 
        k, 
        i2, 
        j2, 
        jj, 
        im, 
        jm,
        *NearNeighb = NULL;
   REAL DMin, 
        x, 
        xx,
        *DissimNearNeighb = NULL,
        *LDDissim = NULL,
        *membr = NULL;
   BOOL *Flag = NULL;

   /* Indicate agglomerable object/clusters                             */
   Flag  = (BOOL *)malloc(NVec * sizeof(BOOL));
   /* Current nearest neighbour storage                                 */
   NearNeighb    = (int  *)malloc(NVec * sizeof(int));
   /* Cluster cardinalities                                             */
   membr = (REAL *)malloc(NVec * sizeof(REAL));
   /* Dissimilarity of nearest neighbour                                */
   DissimNearNeighb = (REAL *)malloc(NVec * sizeof(REAL));
   /* Stores dissimilarities in lower half diagonal                     */
   LDDissim  = (REAL *)malloc(NVec*(NVec-1)/2 * sizeof(REAL));

   /* Check allocations                                                 */
   if(Flag             == NULL || 
      NearNeighb       == NULL || 
      membr            == NULL || 
      DissimNearNeighb == NULL || 
      LDDissim         == NULL)
   {
      if(Flag!=NULL)             free(Flag);
      if(NearNeighb!=NULL)       free(NearNeighb);
      if(membr!=NULL)            free(membr);
      if(DissimNearNeighb!=NULL) free(DissimNearNeighb);
      if(LDDissim!=NULL)         free(LDDissim);

      return(FALSE);
   }
   

   /* For all arrays, move pointer back one so we count FORTRAN-style
      from 1 rather than from 0
   */
   Flag--;                  /* Local arrays                             */
   NearNeighb--;
   membr--;
   DissimNearNeighb--;
   LDDissim--;

   crit--;                  /* Passed parameter arrays                  */
   ib--;
   ia--;
   for(i=0; i<NVec; i++)
      (data[i])--;
   data--;
   

   /* Initializations                                                   */
   for(i=1; i<=NVec; i++) 
   {
      membr[i] = (REAL)1.0;
      Flag[i]  = TRUE;
   }
   NClusters = NVec;
   
   /* Construct dissimilarity matrix                                    */
   for(i=1; i<=NVec-1; i++) 
   {
      for(j=i+1; j<=NVec; j++) 
      {
         ind = IOFFSET(NVec, i, j);
         LDDissim[ind] = (REAL)0.0;
         for(k=1; k<=VecDim; k++) 
         {
            LDDissim[ind] += (data[i][k] - data[j][k]) * 
                             (data[i][k] - data[j][k]);
         }
         
         /* For the case of the min. var. method where merging criteria 
            are defined in terms of variances rather than distances. 
         */
         if (ClusterMethod == 1) 
         {
            LDDissim[ind] /= (REAL)2.0;
         }
      }
   }
   
   /* Carry out an agglomeration - first create list of near neighbours */
   for(i=1; i<=NVec-1; i++) 
   {
      DMin = INF;
      for(j=i+1; j<=NVec; j++) 
      {
         ind = IOFFSET(NVec, i, j);
         if (LDDissim[ind] < DMin) 
         {
            DMin = LDDissim[ind];
            jm = j;
         }
      }
      NearNeighb[i] = jm;
      DissimNearNeighb[i] = DMin;
   }
   
   
   /* Next, determine least dissimilar using list of near neighbours    */
   do
   {
      DMin = INF;
      for(i=1; i<=NVec-1; i++) 
      {
         if(Flag[i] && (DissimNearNeighb[i] < DMin))
         {
            DMin = DissimNearNeighb[i];
            im   = i;
            jm   = NearNeighb[i];
         }
      }
      NClusters--;
      
      /* This allows an agglomeration to be carried out                 */
      i2 = MIN(im,jm);
      j2 = MAX(im,jm);
      ia[NVec   - NClusters] = i2;
      ib[NVec   - NClusters] = j2;
      crit[NVec - NClusters] = DMin;
      
      /* Update dissimilarities from new cluster                        */
      Flag[j2] = FALSE;
      DMin     = INF;
      for(k=1; k<=NVec-1; k++) 
      {
         if(Flag[k] && (k != i2))
         {
            x = membr[i2] + membr[j2] + membr[k];
            
            if (i2 < k) 
               ind1 = IOFFSET(NVec, i2, k);
            else 
               ind1 = IOFFSET(NVec, k, i2);
            
            if (j2 < k)
               ind2 = IOFFSET(NVec, j2, k);
            else
               ind2 = IOFFSET(NVec, k, j2);
            
            ind3 = IOFFSET(NVec, i2, j2);
            xx   = LDDissim[ind3];
            
            switch(ClusterMethod)
            {
            case 1:
               /*  Ward's minimum variance method                       */
               LDDissim[ind1] = (membr[i2] + membr[k]) * LDDissim[ind1] + 
                                (membr[j2] + membr[k]) * LDDissim[ind2] -
                                membr[k] * xx;
               LDDissim[ind1] /= x;
               break;
            case 2:
               /*  Single link method                                   */
               LDDissim[ind1] = MIN(LDDissim[ind1], LDDissim[ind2]);
               break;
            case 3:
               /*  Complete link method                                 */
               LDDissim[ind1] = MAX(LDDissim[ind1], LDDissim[ind2]);
               break;
            case 4:
               /*  Average link (or group average) method               */
               LDDissim[ind1] = (membr[i2] * LDDissim[ind1] + 
                                 membr[j2] * LDDissim[ind2]) / 
                                (membr[i2] + membr[j2]);
               break;
            case 5:
               /*  McQuitty's method                                    */
               LDDissim[ind1] = LDDissim[ind1] * (REAL).5 + 
                                LDDissim[ind2] * (REAL).5;
               break;
            case 6:
               /*  Median (Gower's) method                              */
               LDDissim[ind1] = LDDissim[ind1] * (REAL).5 + 
                                LDDissim[ind2] * (REAL).5 - 
                                xx         * (REAL).25;
               break;
            case 7:
               /*  Centroid method                                      */
               LDDissim[ind1] = (membr[i2] * LDDissim[ind1] + 
                                 membr[j2] * LDDissim[ind2] - 
                                 membr[i2] * membr[j2]  * xx / 
                                 (membr[i2] + membr[j2])) / 
                                (membr[i2] + membr[j2]);
               break;
            }
            
            if((i2 <= k) && (LDDissim[ind1] < DMin))
            {
               DMin = LDDissim[ind1];
               jj = k;
            }
         }
      }
      
      membr[i2] += membr[j2];
      DissimNearNeighb[i2] = DMin;
      NearNeighb[i2] = jj;
      
      /* Update list of nearest neighbours as required.                 */
      for(i=1; i<=NVec-1; i++) 
      {
         if(Flag[i]) 
         {
            if(NearNeighb[i]==i2 || NearNeighb[i]==j2) 
            {
               /* Redetermine nearest neighbour of I                    */
               DMin = INF;
               for(j=i+1; j<=NVec; j++) 
               {
                  ind = IOFFSET(NVec, i, j);
                  if(Flag[j] && (i!=j) && (LDDissim[ind] < DMin))
                  {
                     DMin = LDDissim[ind];
                     jj = j;
                  }
               }
               NearNeighb[i] = jj;
               DissimNearNeighb[i] = DMin;
            }
         }
      }
   }  while(NClusters>1);
   /* Repeat previous steps until N-1 agglomerations carried out.       */

   free(++DissimNearNeighb);
   free(++LDDissim);
   free(++membr);
   free(++NearNeighb);
   free(++Flag);

   data++;
   for(i=0; i<NVec; i++)
      (data[i])++;

   return(TRUE);
}


/************************************************************************/
/*>int **ClusterAssign(FILE *fp, int NVec, int *ia, int *ib, REAL *crit, 
                       int lev, int *iorder, REAL *critval, int *height,
                       LABEL *LabelList)
   ---------------------------------------------------------------------
   Input:   FILE *fp          File for output (or NULL)
            int  NVec         Number of vectors
            int  ia[NVec]     \
                 ib[NVec]     | History of allomerations
            REAL crit[NVec]   /
            int  lev          Number of clusters in largest partition
            LABEL *LabelList  Linked list of labels for the clusters
   Output:  int  iorder[lev]  \
            REAL critval[lev] | vectors describing the dendrogram
            int  height[lev]  /
   Returns: int  **           NVec X lev array describing clusters

   Given a HIERARCHIC CLUSTERING, described as a sequence of    
   agglomerations, derive the assignments into clusters for the 
   top LEV-1 levels of the hierarchy.                           
   Prepare also the required data for representing the          
   dendrogram of this top part of the hierarchy.                

   Pick out the clusters which the N objects belong to, 
   at levels N-2, N-3, ... N-LEV+1 of the hierarchy. 
   The clusters are identified by the lowest seq. no. of 
   their members. 
   There are 2, 3, ... LEV clusters, respectively, for the 
   above levels of the hierarchy. 

   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.     
   22.06.95 Size of hvals array should be [lev+2], not [lev]
   26.06.95 Fixed potential off-bottom of array accesses
   10.06.96 Added LabelList parameter; prints labels
*/
int **ClusterAssign(FILE *fp, int NVec, int *ia, int *ib, REAL *crit, 
                    int lev, int *iorder, REAL *critval, int *height,
                    LABEL *LabelList)
{
   int   ilev, 
         i, 
         j, 
         k, 
         level, 
         icl, 
         NClusters, 
         loc,
         *hvals = NULL,
         **clusters = NULL;
   BOOL  BreakOut;
   LABEL *l;

   /* Allocate memory                                                   */
   clusters = (int **)Array2D(sizeof(int),NVec,lev);
   hvals  = (int *)malloc((2+lev) * sizeof(int));

   /* Check allocations; free and return if failed                      */
   if(clusters==NULL || hvals==NULL)
   {
      if(clusters!=NULL)
         FreeArray2D((char **)clusters, NVec, lev);
      if(hvals!=NULL)
         free(hvals);
      return(NULL);
   }

   /* For all arrays, move pointer back one so we count FORTRAN-style
      from 1 rather than from 0
   */
   crit--;
   ia--;
   ib--;
   height--;
   critval--;
   iorder--;
   hvals--;
   for(i=0; i<NVec; i++)
      (clusters[i])--;
   clusters--;


   hvals[1] = 1;
   hvals[2] = ib[NVec - 1];
   loc = 3;
   for(i=NVec-2; i>=NVec-lev && i>0; i--) 
   {
      for(j=1,BreakOut=FALSE; j<=loc-1; j++) 
      {
         if(ia[i] == hvals[j]) 
         {
            BreakOut=TRUE;
            break;
         }
      }
      if(!BreakOut)
      {
         hvals[loc] = ia[i];
         loc++;
      }
      
      for(j=1,BreakOut=FALSE; j<=loc-1; j++) 
      {
         if(ib[i]==hvals[j]) 
         {
            BreakOut=TRUE;
            break;
         }
      }
      if(!BreakOut)
      {
         hvals[loc] = ib[i];
         loc++;
      }
   }

   for(level=NVec-lev; level<=NVec-2; level++) 
   {
      for(i=1; i<=NVec; i++) 
      {
         icl = i;
         for(ilev=1; ilev<=level; ilev++) 
         {
            if(ib[ilev] == icl) 
            {
               icl = ia[ilev];
            }
         }
         NClusters = NVec-level;
         clusters[i][NClusters-1] = icl;
      }
   }

   for(i=1; i<=NVec; i++) 
   {
      for(j=1; j<=lev-1; j++) 
      {
         for(k=2; k<=lev; k++) 
         {
            if(clusters[i][j] == hvals[k]) 
            {
               clusters[i][j] = k;
               break;
            }
         }
      }
   }

   if(fp!=NULL)
   {
      fprintf(fp,"     SEQ NOS 2CL 3CL 4CL 5CL 6CL 7CL 8CL 9CL\n");
      fprintf(fp,"     ------- --- --- --- --- --- --- --- --- ----\n");
      for(i=1, l=LabelList; i<=NVec; i++)
      {
         fprintf(fp,"%11d",i);
         for(j=1; j<=lev-1; j++)
         {
            fprintf(fp,"%4d",clusters[i][j]);
         }
         if(l!=NULL && l->label!=NULL)
            fprintf(fp,"  %s",l->label);
         
         fprintf(fp,"\n");

         /* Step the label pointer on if there is one                   */
         if(l!=NULL)
            NEXT(l);
      }
   }
   
   /* Determine an ordering of the LEV clusters (at level LEV-1) 
      for later representation of the dendrogram. 
      These are stored in IORDER. 
      Determine the associated ordering of the criterion values 
      for the vertical lines in the dendrogram. 
      The ordinal values of these criterion values may be used in 
      preference, and these are stored in HEIGHT. 
      Finally, note that the LEV clusters are renamed so that they 
      have seq. nos. 1 to LEV.
   */ 
   
   iorder[1] = ia[NVec - 1];
   iorder[2] = ib[NVec - 1];
   critval[1] = (float)0.;
   critval[2] = crit[NVec - 1];
   height[1] = lev;
   height[2] = lev - 1;
   loc = 2;
   for(i=NVec-2; i>=NVec-lev+1; i--) 
   {
      for(j=1; j<=loc; j++) 
      {
         if(ia[i] == iorder[j]) 
         {
            /* Shift rightwards and insert IB(I) beside IORDER(J)       */
            for(k=loc+1; k>=j+1; k--) 
            {
               iorder[k]  = iorder[k-1];
               critval[k] = critval[k-1];
               height[k]  = height[k-1];
            }
            iorder[j + 1]  = ib[i];
            critval[j + 1] = crit[i];
            height[j + 1]  = i - (NVec - lev);
            loc++;
         }
      }
   }

   for(i=1; i<=lev; i++) 
   {
      for(j=1; j<=lev; j++) 
      {
         if(hvals[i]==iorder[j]) 
         {
            iorder[j] = i;
            break;
         }
      }
   }

   /* Fix iorder[] array to give the correct numbers along the bottom   */
   iorder[1] = 1;
   iorder[2] = 2;
   
   for(j=2; j<=lev-1; j++)
   {
      int parent;
      
      for(i=1; i<=NVec; i++)
      {
         if(clusters[i][j] == j+1)
         {
            /* This is the new cluster; see what it's parent was        */
            parent = clusters[i][j-1];

            /* Insert this cluster number to the right of the parent in 
               iorder
            */
            InsertIorder(iorder,lev,j+1,parent);
            
            break;
         }
      }
      
   }
   
   /* Fix pointers into arrays                                          */
   clusters++;
   for(i=0; i<NVec; i++)
      (clusters[i])++;

   /* Free temporary storage                                            */
   free(++hvals);

   /* Return cluster array                                              */
   return(clusters);
}


/************************************************************************/
/*>BOOL InsertIorder(int *iorder, int lev, int cluster, int parent)
   ----------------------------------------------------------------
   Inserts the number cluster into iorder to the right of where parent
   is found. There are lev positions in iorder which is NUMBERED FROM 1

   03.07.95 Original   By: ACRM
*/
BOOL InsertIorder(int *iorder, int lev, int cluster, int parent)
{
   int i,j;
   
   iorder++;
   
   for(i=0; i<lev; i++)
   {
      if(iorder[i]==parent)
      {
         if(i==(lev-1))
            return(FALSE);
         
         for(j=lev-1; j>=i+2; j--)
            iorder[j] = iorder[j-1];
         
         iorder[i+1] = cluster;
         
         return(TRUE);
      }
   }
   return(FALSE);
}

/************************************************************************/
/*>char **ClusterDendogram(FILE *fp, int lev, int *iorder, int *height, 
                           REAL *critval) 
   ----------------------------------------------------------------------
   Input:   FILE *fp          Output file pointer (or NULL)
            int  lev          Number of clustering levels to display
            int  iorder[lev]  Ordering of objects along the bottom of the
                              dendrogram     
            int  height[lev]  Height of the vertical above each object, 
                              in ordinal values   
            REAL critval[lev] Height in real values
   Returns: char **           The dendogram [lev*3][lev*3]
                                                  
   Construct a dendrogram of the top lev levels of a hierarchic 
   clustering.                       

   20.06.95 Original By: ACRM
            Based on FORTRAN code by F. Murtagh, ESA/ESO/STECF, Garching,
            February 1986 available in STATLIB.     
*/
char **ClusterDendogram(FILE *fp, int lev, int *iorder, int *height, 
                        REAL *critval)
{
   int  i,  j,  k,  l,
        i2, j2, i3, 
        ic, idum;
   char **out;
   
   if((out = (char **)Array2D(sizeof(char),lev*3,lev*3))==NULL)
   {
      return(NULL);
   }
   
   /* Blank the dendogram                                               */
   for(i=0; i<lev*3; i++)
   {
      for(j=0; j<lev*3; j++)
      {
         out[i][j]=BLANK;
      }
   }
   
   /* Build the dendogram                                               */
   for(i=3; i<=lev*3; i+=3)
   {
      i2=i/3;
      
      j2=(lev*3+1)-3*height[i2-1];
      for(j=lev*3; j>=j2; j--)
      {
         out[j-1][i-1]=UP;
      }
      
      for(k=i; k>=3; k--)
      {
         i3=(int)((k+2)/3);
         if(((lev*3+1)-height[i3-1]*3)<j2)
            break;
         out[j2-1][k-1]=ACROSS;
      }
   }
   
   /* Print the dendogram                                               */
   if(fp != NULL)
   {
      ic=3;
      for(i=1; i<=lev*3; i++)
      {
         if(i==ic+1)
         {
            idum=ic/3;
            idum=lev-idum;
            for(l=1; l<=lev; l++)
            {
               if(height[l-1]==idum)
                  break;
            }
            idum=l;
            
            fprintf(fp,"         %12.2g    ",critval[idum-1]);
            for(j=0; j<lev*3; j++)
            {
               fprintf(fp,"%c",out[i-1][j]);
            }
            fprintf(fp,"\n");
            
            ic+=3;
         }
         else
         {
            fprintf(fp,"                         ");
            for(j=0; j<lev*3; j++)
            {
               fprintf(fp,"%c",out[i-1][j]);
            }
            fprintf(fp,"\n");
         }
      }
      
      fprintf(fp,"\n                         ");
      for(i=0; i<lev; i++)
         fprintf(fp,"%3d",iorder[i]);
      fprintf(fp,"\n\n");
      
      fprintf(fp,"              CRITERION        CLUSTERS 1 TO LEV\n");
      fprintf(fp,"              VALUES.      (TOP LEV-1 LEVELS OF \
HIERARCHY).\n");
   }

   return(out);
}


/************************************************************************/
/*>BOOL FindMedian(FILE *fp, int **clusters, REAL **data, int NVec, 
                   int VecDim, int NClus, int ClusNum)
   ----------------------------------------------------------------
   For NClus clusters, find the median of cluster ClusNum and then find
   the vector closest to the median. 
   Print this vector.

   04.07.95 Original    By: ACRM
   06.07.95 Fixed return type to BOOL
   02.08.95 Corrected [NClus] tp [NClus-1]
*/
BOOL FindMedian(FILE *fp, int **clusters, REAL **data, int NVec, 
                int VecDim, int NClus, int ClusNum)
{
   int      i, j,
            best;
   REAL     *minval, 
            *maxval,
            *medval,
            mindist,
            dist;
   BOOL     Done = FALSE;

   /* Allocate arrays to store min and max values in each dimension     */
   if((minval=(REAL *)malloc(VecDim*sizeof(REAL)))==NULL)
      return(FALSE);
   if((maxval=(REAL *)malloc(VecDim*sizeof(REAL)))==NULL)
   {
      free(minval);
      return(FALSE);
   }

   /* We just use the same storage space for the median values          */
   medval = minval;

   /* Find the min and max values in each dimension                     */
   for(i=0, Done=FALSE; i<NVec; i++)
   {
      if(clusters[i][NClus-2] == ClusNum)
      {
         /* On first item, just copy in the data                        */
         if(!Done)
         {
            for(j=0; j<VecDim; j++)
               minval[j] = maxval[j] = data[i][j];
            Done = TRUE;
         }
         else
         {
            for(j=0; j<VecDim; j++)
            {
               if(data[i][j] < minval[j])
                  minval[j] = data[i][j];
               if(data[i][j] > maxval[j])
                  maxval[j] = data[i][j];
            }
         }
      }
   }

   /* Now store the median values                                       */
   for(j=0; j<VecDim; j++)
      medval[j] = (minval[j] + maxval[j]) / (REAL)2.0;
      
   /* Now run through again and find which is closest to the medval     */
   for(i=0, Done=FALSE; i<NVec; i++)
   {
      if(clusters[i][NClus-2] == ClusNum)
      {
         if(!Done)
         {
            best    = i;
            mindist = (REAL)0.0;
            for(j=0; j<VecDim; j++)
               mindist += (data[i][j] - medval[j]) *
                          (data[i][j] - medval[j]);
            
            Done = TRUE;
         }
         else
         {
            dist = (REAL)0.0;
            for(j=0; j<VecDim; j++)
               dist += (data[i][j] - medval[j]) *
                       (data[i][j] - medval[j]);
            if(dist < mindist)
            {
               mindist = dist;
               best    = i;
            }
         }
      }
   }

   /* Free up the arrays                                                */
   free(minval);
   free(maxval);

   /* Print the `best' vector                                           */
   for(j=0; j<VecDim; j++)
      fprintf(fp, "%f ",data[best][j]);
   fprintf(fp,"\n");

   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   21.06.95 Original    By: ACRM
   29.06.95 V1.2
   03.07.95 V1.3
   04.07.95 V1.4 Added -n
   06.07.95 V1.5
   28.07.95 V1.6 Added -l
   02.08.95 V1.7  
   11.09.95 V1.8
   10.06.96 V2.0
   20.12.11 V2.1
*/
void Usage(void)
{
   fprintf(stderr,"\nCluster V2.1 (c) 1995 Dr. Andrew C.R. Martin, \
UCL.\n\n");

   fprintf(stderr,"Usage: cluster [-t] [-d] [-l <nclus>] [-n <nclus>] \
[-m <type>] [infile [outfile]]\n");
   fprintf(stderr,"       -t Do not display cluster table\n");
   fprintf(stderr,"       -d Do not display dendogram\n");
   fprintf(stderr,"       -l Limit printing of clusters & dendogram to \
specified number\n");
   fprintf(stderr,"       -n Print near-median vector for specified \
number of clusters\n");
   fprintf(stderr,"       -m Cluster type. This may be a number or word \
as follows:\n");
   fprintf(stderr,"          1 Ward                Ward's minimum \
variance method (default)\n");
   fprintf(stderr,"          2 Single              Single linkage\n");
   fprintf(stderr,"          3 Complete            Complete linkage\n");
   fprintf(stderr,"          4 Average (or Group)  Group average\n");
   fprintf(stderr,"          5 McQuitty            McQuitty's method\n");
   fprintf(stderr,"          6 Median (or Gower)   Gower's median \
method\n");
   fprintf(stderr,"          7 Centroid            Centroid method\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Performs cluster analysis on an arbitrary numerical \
dataset. Each vector\n");
   fprintf(stderr,"to be clustered must appear as a single record in the \
file. No description\n");
   fprintf(stderr,"of the file size is necessary and blank lines are \
ignored. Comments may\n");
   fprintf(stderr,"be placed in the file on lines starting with a ! or \
a #\n\n");
   fprintf(stderr,"As of V2.0, a label may be appended to each vector \
(row) in the data file.\n");
   fprintf(stderr,"This label must consist of a # character followed by \
an arbitrary string.\n");
   fprintf(stderr,"Labels will then be appended to the clustering \
table.\n\n");

   fprintf(stderr,"Clustering code based on code by F. Murtagh, \
ESA/ESO/STECF, Garching.\n");
   fprintf(stderr,"The original FORTRAN code is available from STATLIB \
FTP sites.\n\n");
}


