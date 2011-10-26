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

#define RANGE_LUMINANCE 255 /* Range of image luminance values */
#define PERMS       	0644     /* File acess permits:RW for the users and R for the others */
#define RANGEY      	 255     /* Range level of luminance */ 

#define Clip1(a)            ((a)>5?5:((a)<-5?-5:(a)))

int* calcula_media(int **set_treino);
void ajuda( char *prgname );
void read_header_pgm (int *ysize, int *xsize, char *Fname);
void read_f_pgm (int **pelimg, int *ysize, int *xsize, char *Fname);
int **intmatrix(int nr, int nc);
float **floatmatrix(int nr, int nc);
void cria_dicionario(void);
float quad_err(int ind_codebook, int ind_set_treino);
void escreve_codebook(void);

int **Image;
int **Image_out;
int ysize[1], xsize[1];           /* The dimentions of the original image */

int BsizeX;
int BsizeY;

int num_vectores;
float num_vectores_total;
int **set_treino;

int *vector_medias;
int *vector_perturbacao;

int num_codewords;
int num_final;
int **codebook;
int *cont_final;
float **medias_final;
int THRESHOLD;

FILE *pointf_out;
char *inname, *outname;			//Nome dos ficheiros de input e output

int info=0;

//******************************************************************************
//*                                                                            *
//*     Funcao main                                                            *
//*                                                                            *
//*                                                                            *
//******************************************************************************
int main(int argc, char *argv[])
{

	int i, j, i1, j1, n;
	int contador;
	float media=0;

	clock_t start, end;
	double elapsed=0;
	start = clock();

#ifdef VARS_NOT_USED
	int aux2;
	int aux;
#endif /* VARS_NOT_USED */
 	
	THRESHOLD=10;				//Criterio de paragem para convergencia do algoritmo LBG

	srand48( (unsigned int)time( NULL ) );

    	//Verifica argumentos
   	if( argc < 6 ){
        	ajuda(argv[0]);
        return 0;
    	}

	num_final=atoi(argv[1]);

	BsizeY = atoi(argv[2]);
	BsizeX = atoi(argv[3]);

	outname = argv[4];

	info=atoi(argv[5]);

	printf("\n--------------------------------------------------------");
	printf("\n Criando dicionario: %s", outname);
	printf("\n Numero de codewords a gerar: %d", num_final);
	printf("\n Blocos %dx%d", BsizeY, BsizeX);


	num_vectores_total=0;
	
	//Le as varias imagens de treino
	for(n=6; n<argc; n++){					
		inname  = argv[n];

		printf("\n--------------------------------------------------------");
		printf("\n Lendo imagem %s", inname);
		read_header_pgm (ysize, xsize, inname);     /* Reads the PGM file and returns the picture size */

		free(Image);
  		Image=intmatrix(*ysize, *xsize);

  		read_f_pgm (Image, ysize, xsize, inname);    /* Reads the PGM file and stores the image in pely */

		printf("\n Tamanho (%dx%d) %d", *xsize, *ysize, ((*xsize)*(*ysize)));
	
		num_vectores=(*ysize/BsizeY)*(*xsize/BsizeX);
  		printf("\n Numero de blocos %dx%d: %d", BsizeY, BsizeX, num_vectores);

		media=0;
		//Calcula a media das luminancias da imagem
		for(i=0; i<*ysize; i++){
    			for(j=0; j<*xsize; j++){
				media+=Image[i][j];
			}
		}
		media=media/((*xsize)*(*ysize));

		
		printf("\n Media: %.0f", media);

		//Subtrai a media a todos os pixels
		for(i=0; i<*ysize; i++){
    			for(j=0; j<*xsize; j++){
				Image[i][j]-=media;
			}
		}
		//-----------------------------------------

		
		//Aloca um vector para receber todos os elementos de treino  
  		set_treino=(int **)realloc(set_treino,(int)(num_vectores_total+num_vectores)*sizeof(int *));
  		if (!set_treino){
    		printf("intmatrix() - allocation failure 1 \n");
    		exit(1);
  		}
  
  		for(i=num_vectores_total;i<num_vectores_total+num_vectores;i++) {
    			set_treino[i]=(int *)malloc((unsigned) (BsizeY*BsizeX)*sizeof(int));
    			if (!set_treino[i]){
      				printf("intmatrix() - allocation failure 2 \n");
      				exit(1);
    			}      
  		}
		//-----------------------------------------

		//Lê os vectores de treino
  		for(i=0, contador=num_vectores_total; i<*ysize; i+=BsizeY)
    			for(j=0; j<*xsize; j+=BsizeX, contador++){
				for(i1=0; i1<BsizeY; i1++)
	    				for(j1=0; j1<BsizeX; j1++){
						set_treino[contador][j1+(i1*BsizeX)]=Image[i+i1][j+j1];
	      				}
	 		}

		num_vectores_total+=num_vectores;
		//-----------------------------------------
	}


	
	printf("\n--------------------------------------------------------");
	printf("\n Numero total de blocos lidos : %.0f", num_vectores_total);
	printf("\n Processando                  :     ");
	fflush(stdout);
	//-----------------------------------------
	num_vectores=num_vectores_total;

	if(num_vectores_total<num_final){
		printf("\n--------------------------------------------------------");
		printf("\n O numero de vectores a gerar e superior ao numero de vectores da imagem");
		printf("\n Introduza mais imagens ou reduza o numero de vectores a gerar");
		printf("\n--------------------------------------------------------\n\n");
		exit(1);
	}


	
	cria_dicionario();

	if(info==1){
	  for(contador=0; contador<num_final; contador++){
		printf("\nVector[%d]: %d elementos", contador, cont_final[contador]);
	  }

	  printf("\n");
	}


	end = clock();
 	elapsed += ((double) (end - start)) / CLOCKS_PER_SEC;
 	start = clock();

	printf("\n--------------------------------------------------------");
	printf("\n Construccao do dicionario concluida");
	printf("\n Tempo total de execucao      : %9.3f segundos", elapsed);
	printf("\n--------------------------------------------------------\n\n");

	escreve_codebook();

  	return EXIT_SUCCESS;
}
//Fim da funcao main








void cria_dicionario(void){
	int i, j, n;
	int contador;
	int indice;
	int limite;
	float distorcao;
	float aux;
	float distorcao_total;
	float distorcao_anterior=FLT_MAX;

	int new_num;
	int max_elementos;
	int elemento_a_dividir;
	int elementos_nulos=0;
	int *dividido;

#ifdef VARS_NOT_USED
	int interacoes;
	int auxaux=0;
	int percentagem=0;
#endif /* VARS_NOT_USED */

	//------------------------------------------------------
	codebook=intmatrix(num_final, BsizeY*BsizeX);

	medias_final=floatmatrix(num_vectores, BsizeY*BsizeX);

	cont_final=(int *)calloc(num_final,sizeof(int));
  	if (!cont_final){
    		printf("floatmatrix() - allocation failure 1 \n");
    		exit(1);
  	}
	//------------------------------------------------------
	dividido=(int *)calloc(num_final,sizeof(int));
  	if (!cont_final){
    		printf("floatmatrix() - allocation failure 1 \n");
    		exit(1);
  	}

	codebook[0]=calcula_media(set_treino);	
	num_codewords=1;


	do{
	  printf("\b\b\b\b%3d%%", ((num_codewords*100)/num_final));
	  fflush(stdout);
	
	  for(i=0; i<num_codewords; i++){
		dividido[i]=0;
	  }	

	  new_num=0;

	  for(n=0; n<num_codewords; n++){
		for(i=0; i<BsizeX*BsizeY; i++){
			limite=ceil((255/(1+(sqrt(num_codewords))))/2);
			codebook[n+num_codewords][i]=codebook[n][i]+mrand48()%(limite); //limita-se amplitude do vector de perturbacao a medida que o numero de outputs aumenta
			if(codebook[n+num_codewords][i]>255){
				codebook[n+num_codewords][i]=255;
			}else{
				if(codebook[n+num_codewords][i]<-255){
					codebook[n+num_codewords][i]=-255;
				}
			}	
 		}
		new_num++;
		if(new_num+num_codewords>=num_final){
			break;
		}

	  }

	  num_codewords+=new_num;

	  if(info==1){
	    printf("\n %d\n", num_codewords);
	    fflush(stdout);
	  }


	  //----------------------------------------------------------
	  for(;;){
	    do{

		for(i=0; i<num_codewords; i++){	
			cont_final[i]=0;
			for(j=0; j<BsizeX*BsizeY; j++){
				medias_final[i][j]=0;
			}		
		}


	    distorcao_total=0;
	    indice=0;

	    for(contador=0; contador<num_vectores; contador++){		//Para todos os vectores do set de treino
	      distorcao=FLT_MAX;
		for(i=0; i<num_codewords; i++){				//Varre todos os elementos do codebook
			aux=quad_err(i,contador);
			//printf("Distorcao[%d]= %f\t", i, aux);
			if(aux<distorcao){
				distorcao=aux;
				indice=i;
			}			
		}
		//printf("Indice %d\t", indice);							//Quando termina este ciclo, é determinado o codeword de distorcao minima
		distorcao_total+=distorcao;
		//printf("Distorcao: %f, Distorcao_total: %f\n",distorcao, distorcao_total);
		cont_final[indice]++;		
		for(j=0; j<BsizeX*BsizeY; j++){
			medias_final[indice][j]+=set_treino[contador][j];
		}		
	  }


	  for(i=0; i<num_codewords; i++){
		for(j=0; j<BsizeX*BsizeY; j++){
			if(cont_final[i]!=0){
				medias_final[i][j]=medias_final[i][j]/cont_final[i];
				codebook[i][j]=(int)medias_final[i][j];
			}
		}
	  }
		
	  //Procura por elementos sem vectores associados
	  elementos_nulos=0;
	  for(i=0; i<num_codewords; i++){
		elemento_a_dividir=-1;
		if(cont_final[i]==0){

			if(info==1){
			  printf("\nElemento vazio [%d]", i);
			}
			elementos_nulos++;

			//Procura o elemento com mais vectores associados
			max_elementos=0;
			for(j=0; j<num_codewords; j++){
				if((cont_final[j]>max_elementos)&&(dividido[j]<4)){		//Salvaguarda situacoes em que a divisao do vector nao traz alteracoes, para evitar loops
					max_elementos=cont_final[j];
					elemento_a_dividir=j;
				}
			}

			if(elemento_a_dividir>=0){
			  dividido[elemento_a_dividir]++;						//Indica que o vector ja foi dividido uma vez
			  if(info==1){
				printf("\telemento[%d] = elemento[%d]", i, elemento_a_dividir);
			  }
			  for(n=0; n<BsizeX*BsizeY; n++){
				codebook[i][n]=codebook[elemento_a_dividir][n]+mrand48()%5;
			  }
			}else{
				elementos_nulos=0;		//Ja nao ha ninguem para dividir
			}
		}
	}

	//-------------------------------------------------
	}while(elementos_nulos>0);//Nao existe nenhum elemento sem vectores associados

	if(info==1){
		printf("\nDistorcao_total: %f", distorcao_total);
	}

	//Criterio de paragem
	if((distorcao_anterior-distorcao_total)<THRESHOLD){
		break;
	}
	distorcao_anterior=distorcao_total;

	}

    }while(num_codewords<num_final);
    printf("\n Distorcao total final        : %.0f", distorcao_total);
}




//******************************************************************************
//*                                                                            *
//*     Escreve o codebook gerado para o ficheiro especificado                 *
//*                                                                            *
//******************************************************************************
void escreve_codebook(void){
int i, j;

	pointf_out=fopen(outname, "w");

  	if(pointf_out==NULL)
    	{
      		fprintf(stderr,"Impossivel abrir ficheiro: %s\n\n",outname);
		exit(1);
    	}

	fprintf(pointf_out, "%d\n", num_final);
	fprintf(pointf_out, "%d\n", BsizeX);
	fprintf(pointf_out, "%d\n", BsizeY);
	

	for(i=0; i<num_codewords; i++){
		for(j=0; j<BsizeX*BsizeY; j++){
			fprintf(pointf_out, "%d\t", codebook[i][j]);
		}
		fprintf(pointf_out, "\n");
	}
fclose(pointf_out);
}




//******************************************************************************
//*                                                                            *
//*     Calcula o erro quadrático entre um vectro do set de treino             *
//*	e um vector do codebook                                                *
//*                                                                            *
//******************************************************************************
float quad_err(int ind_codebook, int ind_set_treino){
	int i;
	float tmp=0;

	for(i=0; i<BsizeX*BsizeY; i++){
		tmp+=((codebook[ind_codebook][i]-set_treino[ind_set_treino][i])*(codebook[ind_codebook][i]-set_treino[ind_set_treino][i]));
	}

	return tmp;
} /* Enf of quad_err */





//******************************************************************************
//*                                                                            *
//*     Calcula o vector media de um conjunto de vectores                      *
//*                                                                            *
//******************************************************************************
int* calcula_media(int **set_treino){
	int *vector_medias;
	float *acumulador;
	int i, j;

	vector_medias=(int *)calloc(BsizeY*BsizeX,sizeof(int));
	acumulador=(float *)calloc(BsizeY*BsizeX,sizeof(float));

	for(i=0; i<num_vectores; i++){
		for(j=0; j<BsizeX*BsizeY; j++){
			acumulador[j]+=set_treino[i][j];
		}
	}
	
	for(j=0; j<BsizeX*BsizeY; j++){
		vector_medias[j]=acumulador[j]/num_vectores;
	}
	free(acumulador);
	return vector_medias;
}



//******************************************************************************
//*                                                                            *
//*     Apresenta informacao e sintaxe do programa			       *
//*                                                                            *
//******************************************************************************
void ajuda( char *prgname ){
	printf("\n-------------------------------------------------------------------------------------------------------------------------------------------\n");
	printf(" Programa de treino do dicionario baseado no algoritmo Linde-Buzo-Gray com inicializacao por vector medio e vectores perturbacao\n");
	printf("-------------------------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Sintaxe: %s <N� vectores a gerar> <Tamanho X> <Tamanho Y> <Ficheiro de saida> <info> <Imagem de treino 1> <Imagem de treino2> ...\n", prgname);
	printf("-------------------------------------------------------------------------------------------------------------------------------------------\n");
	printf(" info 0-1\n");
	printf("-------------------------------------------------------------------------------------------------------------------------------------------\n\n");
}


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
  
#ifdef VARS_NOT_USED
  char pelimg[1][20];
#endif /*VARS_NOT_USED */

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




/*
 * read_f_pgm: READS A PGM FILE
 * usage: read_f_pgm (unsigned char **pelimg, int *ysize,
 * 				 		int *xsize, char *Fname) 
 *
 * Inputs:
 * Fname - File name
 * Returns 
 * - the image in char **pellimg                                              
 * - number of rows 
 * - number of columns
 */
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



/*************************************************************************************/
/*                                                                                   */
/*  FLOATMATRIX - Allocates memory for a matrix of variables of type float           */
/*                                                                                   */
/*  Inputs:                                                                          */
/*    number of rows / nunmber of columnc                                            */
/*  Returns a poiter to a int matrix (float **)                                        */
/*                                                                                   */
/*************************************************************************************/
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


