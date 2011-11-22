/***************************************************************************
 *   Copyright (C) 2008 by Nelson Carreira Francisco                       *
 *   eng.nelsito@gmail.com                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/***************************************************************************
*   Implementacaoo de um codificadror de imagens baseado em                *
*   Quantificacao vectorial                                                *
*   Nelson Carreira Francisco                                              *
****************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#define RANGE_LUMINANCE 255 /* Range of image luminance values */
#define PERMS       	0644     /* File acess permits:RW for the users and R for the others */
#define RANGEY      	 255     /* Range level of luminance */ 

#define Clip1(a)            ((a)>255?255:((a)<0?0:(a)))


void ajuda( char *prgname );
int **intmatrix(int nr, int nc);
float **floatmatrix(int nr, int nc);
float quad_err(int indice_dic);
void carrega_dicionario(void);
void write_f_pgm ( unsigned char **im_matrix, int nline, int npixel, char *filename);
unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch);
int le_indice(void);
double PSNR (int **origblk, int **cmpblk, int nline, int npixel );
double MSE (int **origblk, int **cmpblk, int nline, int npixel );
void read_header_pgm (int *ysize, int *xsize, char *Fname);
void read_f_pgm (int **pelimg, int *ysize, int *xsize, char *Fname);

/* The BIT BUFFER */
static int buffer;		/* Bits waiting to be input 		*/
static int bits_to_go;		/* Number of bits still in buffer	*/
static int garbage_bits;	/* Number of bits past end-of-file 	*/
long Bits_Count;

void start_inputing_bits();
int input_bit(FILE* file);

int **Image_out;
int ysize[1], xsize[1];           /* The dimentions of the original image */

int BsizeX;
int BsizeY;
int num_codewords;
int bits_indice;
int num_blocos;
int **dic;

int *bloco_original;


FILE *pointf_in;
FILE *pointf_dic;
char *inname, *outname, *dic_name;			//Nome dos ficheiros de input e output
char *inimage;
int **Image;
//******************************************************************************
//*                                                                            *
//*     Funcao main                                                            *
//*                                                                            *
//*                                                                            *
//******************************************************************************
int main(int argc, char *argv[])
{

	int i, j, i1, j1;
	int indice;
	int media;
#ifdef UNUSED_VARS
	int n;
	float distorcao;
	float aux;
	int contador;
#endif /* UNUSED_VARS */

	clock_t start, end;
	double elapsed=0;
	start = clock();
	double psnr, mse;

    	//Verifica argumentos
   	if( argc < 3 ){
        	ajuda(argv[0]);
        return 0;
    	}

        inname  = argv[1];
        dic_name = argv[2];
	outname = argv[3];


	//Carrega dicionario
	carrega_dicionario();
	bits_indice=ceil(log(num_codewords)/log(2));

	pointf_in=fopen(inname, "r");
  	if(pointf_in==NULL)
    	{
      		fprintf(stderr,"Impossivel ficheiro: %s\n\n",inname);
		exit(1);
    	}


	printf("\n Ficheiro a descodificar       : %s", inname); 
	fscanf(pointf_in, "%d\n", xsize);
	fscanf(pointf_in, "%d\n", ysize);
	fscanf(pointf_in, "%d\n", &media);

	printf("\n Tamanho (%dx%d)             : %d pixels", *xsize, *ysize, ((*xsize)*(*ysize)));

	Image_out=intmatrix(*ysize, *xsize);

	num_blocos=(*ysize/BsizeY)*(*xsize/BsizeX);
  	printf("\n Total de blocos %dx%d na imagem : %d blocos", BsizeY, BsizeX, num_blocos);
	printf("\n-----------------------------------------------------");

	start_inputing_bits();

 	for(i=0; i<*ysize; i+=BsizeY){
     		for(j=0; j<*xsize; j+=BsizeX){
			indice=le_indice();
 			for(i1=0; i1<BsizeY; i1++){
 	    			for(j1=0; j1<BsizeX; j1++){
 					Image_out[i+i1][j+j1]=dic[indice][i1*BsizeX+j1];
 	      			}
 	  		}
 		}
 	}




	//Soma a média a todos os pixels
	unsigned char **Image_aux=ucmatrix(0, *ysize, 0, *xsize);
	for(i=0; i<*ysize; i++){
    		for(j=0; j<*xsize; j++){
			Image_out[i][j]+=media;
			Image_aux[i][j]=(unsigned char)(Image_out[i][j]);
		}
	}
	//-----------------------------------------

	end = clock();
 	elapsed += ((double) (end - start)) / CLOCKS_PER_SEC;
 	start = clock();

 	write_f_pgm(Image_aux, *ysize, *xsize, outname);
	
	
	printf("\n Imagem descodificada para     : %s", outname);
	printf("\n Tempo total de execucao       : %9.3f segundos", elapsed);
	printf("\n-----------------------------------------------------");


	//Opcionalmente, calcula o PSNR e a MSE se o ficheiro original for fornecido
	if( argc > 4 ){
		inimage  = argv[4];
		read_header_pgm (ysize, xsize, inimage);
		Image=intmatrix(*ysize, *xsize);
		read_f_pgm (Image, ysize, xsize, inimage);
		psnr=PSNR(Image, Image_out, *ysize, *xsize);
		mse=MSE(Image, Image_out, *ysize, *xsize);
		printf("\n Abrindo imagem original       : %s", inimage);
		printf("\n PSNR                          : %f dB", psnr); 
		printf("\n MSE                           : %f", mse);
		printf("\n-----------------------------------------------------");
    	}

  fclose(pointf_in);
  printf("\n\n");
  return EXIT_SUCCESS;
}
//Fim da funcao main



int le_indice(void){
int k;
int bit_readed;
int indice=0;

	for(k=0;k<bits_indice; k++){
		indice=indice<<1;
		bit_readed=input_bit(pointf_in);
		indice=indice|bit_readed;
	}

	return indice;
}






//******************************************************************************
//*                                                                            *
//*     Calcula o erro quadrático entre um vectro do set de treino             *
//*	e um vector do codebook                                                *
//*                                                                            *
//******************************************************************************
void carrega_dicionario(void){
	int i, j;

	pointf_dic=fopen(dic_name, "r");
  	if(pointf_dic==NULL)
    	{
      		fprintf(stderr,"Impossivel abrir dicionario: %s\n\n",dic_name);
		exit(1);
    	}
 
	fscanf(pointf_dic, "%d\n", &num_codewords);
	fscanf(pointf_dic, "%d\n", &BsizeX);
	fscanf(pointf_dic, "%d\n", &BsizeY);

	printf("\n-----------------------------------------------------");
	printf("\n Carregou dicionario %s", dic_name); 
	printf("\n %d blocos de %dx%d pixels", num_codewords, BsizeY, BsizeX);
	printf("\n-----------------------------------------------------");
	fflush(stdout);

	dic=intmatrix(num_codewords, BsizeY*BsizeX);

	for(i=0; i<num_codewords; i++){
		for(j=0; j<BsizeX*BsizeY; j++){
			fscanf(pointf_dic, "%d\t", &dic[i][j]);
		}
		fscanf(pointf_dic, "\n");
	}

	fclose(pointf_dic);
}





//******************************************************************************
//*                                                                            *
//*     Apresenta informacao e sintaxe do programa			       *
//*                                                                            *
//******************************************************************************
void ajuda( char *prgname ){
	printf("\n--------------------------------------------------------------------------------------------------------------\n");
	printf(" Programa de descodificação de imagens baseado em quantização vectorial\n");
	printf("--------------------------------------------------------------------------------------------------------------\n");
    	printf(" Sintaxe: %s <ficheiro origem> <ficheiro dicionario> <ficheiro destino> <imagem original (opcional)>\n\n", prgname);
	printf(" Imagem original permite calcular o PSNR e o MSE para a imagem reconstruida\n");
	printf("--------------------------------------------------------------------------------------------------------------\n\n");
}



/*************************************************************************************/
/*                                                                                   */
/*  INTMATRIX - Allocates memory for a matrix of variables of type int               */
/*                                                                                   */
/*  Inputs:                                                                          */
/*    number of rows / nunmber of columnc                                            */
/*  Returns a poiter to a int matrix (int **)                                        */
/*                                                                                   */
/*************************************************************************************/
int **intmatrix(int nr, int nc)
{
  int i;
  int **m;
  
  m=(int **)malloc((unsigned) (nr)*sizeof(int *));
  if (!m){
    printf("intmatrix() - allocation failure 1 \n");
    exit(1);
  }
  
  for(i=0;i<nr;i++) {
    m[i]=(int *)malloc((unsigned) (nc)*sizeof(int));
    if (!m[i]){
      printf("intmatrix() - allocation failure 2 \n");
      exit(1);
    }      
  }

  return m;
}


float **floatmatrix(int nr, int nc)
{
  int i;
  float **m;
  
  m=(float **)malloc(nr*nc*sizeof(float *));
  if (!m){
    printf("floatmatrix() - allocation failure 1 \n");
    exit(1);
  }
  
  for(i=0;i<nr;i++) {
    m[i]=(float *)calloc(nc, sizeof(float));
    if (!m[i]){
      printf("floatmatrix() - allocation failure 2 \n");
      exit(1);
    }      
  }

  return m;
}




/************************************************************************************/
          /* WRITE OUTPUT LUMINANCE FILE - PGM format*/
/************************************************************************************/
void write_f_pgm ( unsigned char **im_matrix, int nline, int npixel, char *filename)
{
 int i;
 int pointfo;
 char header_pgm[20];
 int npixel_orig, nline_orig;

 npixel_orig=npixel;
 nline_orig=nline;

 if ((pointfo=creat(filename,PERMS))==-1)	/* creates file */
 {
	fprintf(stderr,"Write_f_PGM - cannot create PGM file -> %s",filename);
	exit(1);
 }
 /* Defines the header of file - type pgm */
 *(header_pgm+0)= 'P';
 *(header_pgm+1)= '5';
 *(header_pgm+2)= '\n';
 i=3;
 if(npixel>=1000)
 {
	*(header_pgm+i)= (unsigned char) (npixel/1000+48); 
     	npixel=npixel%1000;
     	i++;
 }
 *(header_pgm+i)= (unsigned char) (npixel/100+48);   i++;
 *(header_pgm+i)= (unsigned char)(npixel%100)/10+48;   i++;
 *(header_pgm+i)= (unsigned char) ((npixel%100)%10+48);   i++;
 *(header_pgm+i)= ' ';  i++;
 if(nline>=1000)
 {
	*(header_pgm+i)= (unsigned char) (nline/1000+48); 
	nline=nline%1000;
	i++;
 }
 *(header_pgm+i)= (unsigned char) (nline/100+48);   i++;
 *(header_pgm+i)= (unsigned char) ((nline%100)/10+48);   i++;
 *(header_pgm+i)= (unsigned char) ((nline%100)%10+48);   i++;
 *(header_pgm+i)='\n';   i++;
 *(header_pgm+i)= (unsigned char) (RANGE_LUMINANCE/100+48);   i++;
 *(header_pgm+i)= (unsigned char) ((RANGE_LUMINANCE%100)/10+48);   i++;
 *(header_pgm+i)= (unsigned char) ((RANGE_LUMINANCE%100)%10+48);   i++;
 *(header_pgm+i)='\n';    i++;  
  write(pointfo, (char *) header_pgm,i);
  for (i=0;i<(nline_orig);i++)
    	write(pointfo, (char *) im_matrix[i],npixel_orig);
  close(pointfo);				/* closes file */
}/*END of write_f_pgm */

/*************************************************************************************/
/*                                                                                   */
/*  UCMATRIX - Allocates memory for a matrix of variables of type unsigned char      */
/*                                                                                   */
/*  Inputs:                                                                          */
/*    1st and last row / 1st and last column                                         */
/*  Returns a poiter to a unsigned char matrix (unsigned char **)                    */
/*                                                                                   */
/*************************************************************************************/
unsigned char **ucmatrix(int nrl,int nrh,int ncl,int nch)
{
  int i;
  unsigned char **m;
  
  m=(unsigned char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
  if (!m){
    printf("ucmatrix() - allocation failure 1 \n");
    exit(1);
  }
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++) {
    m[i]=(unsigned char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
    if (!m[i]){
      printf("ucmatrix() - allocation failure 2 \n");
      exit(1);
    }      
    m[i] -= ncl;
  }
  return m;
}


/* INITIALISE BIT INPUT */
void start_inputing_bits()
{ 
 bits_to_go = 0;						/* Buffer starts out with		*/
 garbage_bits = 0;						/* no bits in it.			*/
 Bits_Count=0;
}


/* INPUT A BIT  */
int input_bit(FILE* file)
{ 
  int t;
  if (bits_to_go==0) 		            		/* read the next byte if no	*/
    {  
      buffer = getc(file);	            		/* bits are left in buffer	*/

      if (buffer==EOF) 
	{
	  garbage_bits +=1;						/* Return arbitrary bits 	*/
        }
      bits_to_go = 8;
    }
  t = buffer&1;									/* Return the next bit from	*/
  buffer >>= 1;									/* the bottom of the byte 	*/
  bits_to_go -= 1;
  Bits_Count ++;
  return(t);
}


/************************************************************************************/
/* Peak Signal Noise Ratio                                                          */
/************************************************************************************/
double PSNR (int **origblk, int **cmpblk, int nline, int npixel )
{
 int i, j;
 double psnr;

  psnr=0.;
  for(j=0;j<nline;j++)
     for(i=0;i<npixel;i++)
         psnr+=((double)*(origblk[j]+i)-(double)*(cmpblk[j]+i)) * ((double)*(origblk[j]+i)-(double)*(cmpblk[j]+i));

   if (psnr==0.0) psnr = 0.000000000000000000000000000000000000000000000000000001 ;
   psnr=10.0*log10((double)(RANGEY*RANGEY)/(psnr/(double)(nline*npixel)));

 return psnr;
}
/* End of Psnr function */


/************************************************************************************/
/* Mean Squared Error                                                               */
/************************************************************************************/
double MSE (int **origblk, int **cmpblk, int nline, int npixel )
{
  int i, j;
  long cnt=0;
  double mse;
  
  mse=0.;
  for(j=0;j<nline;j++)
    for(i=0;i<npixel;i++)
      {
	mse+=((double)*(origblk[j]+i)-(double)*(cmpblk[j]+i)) * ((double)*(origblk[j]+i)-(double)*(cmpblk[j]+i));
	cnt++;
      }
  
  return(mse/cnt);
}
/* End of MSE function */



/************************************************************************************/
/*                                                                                  */
/* READS THE INFORMATION OF A PGM FILE HEADER                                       */
/* usage: read_header_pgm (int *ximg,int *yimg,char *originalfilename);             */
/*                                                                                  */
/************************************************************************************/
void read_header_pgm (int *ysize, int *xsize, char *Fname)
{
  FILE *pointf;
  char dummy[12];
  char aux;
  int lido=0;
  
  pointf=fopen(Fname, "r");
  if(pointf==NULL)
    {
      fprintf(stderr,"\nCannot open input file: %s\n",Fname);
      exit(1);
    }

  fscanf(pointf, "%s\n", dummy);


//*******************************
do{
  aux=getc(pointf);
  if(aux=='#'){
	do{
        aux=getc(pointf);
	}while((aux!='\n')&&(aux!=EOF));

       
  }else{
       ungetc(aux,pointf);
       lido=1;
  }
}while(lido==0);	

        fscanf(pointf," %d %d", xsize, ysize);
//********************************

  fclose(pointf); 		/* closes file */
}

/* End of read_header_pgm function */




/************************************************************************************/
/*                                                                                  */
/* read_f_pgm: READS A PGM FILE                                                     */
/* usage: read_f_pgm (unsigned char **pelimg, int *ysize, int *xsize, char *Fname)  */
/*                                                                                  */
/* Inputs:                                                                          */
/* Fname - File name                                                                */   
/*                                                                                  */
/* Returns 
   - the image in char **pellimg                                              
   - number of rows 
   - number of columns                                                              */
/*                                                                                  */
/************************************************************************************/

void read_f_pgm (int **pelimg, int *ysize, int *xsize, char *Fname)
{
  int i, j;
  FILE *pointf;
  char dummy[15];
  int aux[1];
  char aux1;
  int lido=0;
  
  
  pointf=fopen(Fname, "r");
  if(pointf==NULL)
    {
      fprintf(stderr,"\nCannot open input file: %s\n",Fname);
      exit(1);
    }
  
  fscanf(pointf, "%s\n", dummy);

//*******************************
do{
  aux1=getc(pointf);
  //fflush(stdout);
  if(aux1=='#'){
	do{
        aux1=getc(pointf);
	}while((aux1!='\n')&&(aux1!=EOF));
  }else{
       ungetc(aux1,pointf);
       lido=1;
  }	
}while(lido==0);



        fscanf(pointf," %d %d", xsize, ysize);
//********************************

  fscanf(pointf, "%d\n", aux);

  for (i=0;i<*ysize;i++){
    for (j=0;j<*xsize;j++){
      pelimg[i][j]=(unsigned char) fgetc(pointf);
    }
  }

  fclose(pointf); 		/* closes file */
}
/* End of read_f_pgm function */
