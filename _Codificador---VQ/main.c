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
 *   Implementacaoo de um codificador de imagens baseado em                *
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


void ajuda(char *prgname);
void read_header_pgm(int *ysize, int *xsize, char *Fname);
void read_f_pgm(int **pelimg, int *ysize, int *xsize, char *Fname);
int **intmatrix(int nr, int nc);
float **floatmatrix(int nr, int nc);
float quad_err(int indice_dic);
void carrega_dicionario(void);
double PSNR(int **origblk, int **cmpblk, int nline, int npixel);
double MSE(int **origblk, int **cmpblk, int nline, int npixel);
void escreve_indice(int indice);

void write_f_pgm(int **im_matrix, int nline, int npixel, char *filename);
unsigned char **ucmatrix(int nrl, int nrh, int ncl, int nch);


//Funcoes de escrita de bits para o ficheiro
void output_bit(int bit, FILE* output_file);
void done_outputing_bits(FILE* output_file);

/* THE BIT BUFFER */
static int buffer; /* Bits buffered for output 		*/
static int bits_to_go; /* Number of bits still in buffer	        */
long Bits_Count;


/*
int **Image_orig;
int **Image_out;
int ysize[1], xsize[1];           // The dimensions of the original image 
 */

/*
 * The following variables are loaded from the dictionary  file
 **/
int G_BsizeX;
int G_BsizeY;
int G_num_codewords;


int G_bits_indice;
int **G_dic;

int *bloco_original;

//FILE *pointf_in;		//Ponteiro para o ficheiro (imagem) original
//FILE *pointf_out;		//Ponteiro para o ficheiro de saida

FILE *pointf_in;
FILE *pointf_dic;
FILE *pointf_out;
char *inname, *outname, *dic_name; //Nome dos ficheiros de input e output

//******************************************************************************
//*                                                                            *
//*     Funcao main                                                            *
//*                                                                            *
//*                                                                            *
//******************************************************************************

int main(int argc, char *argv[]) {

    int **Image_orig;
    int **Image_out;
    int ysize[1], xsize[1]; /* The dimensions of the original image */

    int **Image;
    int num_blocos;

    int i, j, i1, j1, n;
    float aux;
    int indice;
#ifdef VARS_NOT_USED
    int contador;
#endif /* VARS_NOT_USED */
    int media;
    float distorcao;
    float acumulador = 0;
    double psnr, mse;
    float taxa;

    clock_t start, end;
    double elapsed = 0;
    start = clock();

    /* Dummy value -- indice is set before its usage */
    indice = -1;

    //Verifica argumentos
    if (argc < 3) {
        ajuda(argv[0]);
        return 0;
    }

    inname = argv[1];
    dic_name = argv[2];
    outname = argv[3];


    //Carrega dicionario
    carrega_dicionario();
    G_bits_indice = ceil(log(G_num_codewords) / log(2));

    bloco_original = (int *) calloc(G_BsizeX*G_BsizeY, sizeof (int));
    if (!bloco_original) {
        printf("intmatrix() - allocation failure 1 \n");
        exit(1);
    }

    //Le imagem a comprimir
    printf("\n Imagem a comprimir            : %s", inname);
    read_header_pgm(ysize, xsize, inname); /* Reads the PGM file and returns the picture size */
    Image_orig = intmatrix(*ysize, *xsize);
    Image = intmatrix(*ysize, *xsize);
    Image_out = intmatrix(*ysize, *xsize);
    read_f_pgm(Image, ysize, xsize, inname); /* Reads the PGM file and stores the image in pely */
    printf("\n Tamanho (%dx%d)             : %d pixels", *xsize, *ysize, ((*xsize)*(*ysize)));

    for (i = 0; i<*ysize; i++) {
        for (j = 0; j<*xsize; j++) {
            Image_orig[i][j] = Image[i][j];
        }
    }

    num_blocos = (*ysize / G_BsizeY)*(*xsize / G_BsizeX);
    printf("\n Total de blocos %dx%d na imagem : %d blocos",
            G_BsizeY, G_BsizeX, num_blocos);
    printf("\n-----------------------------------------------------");

    //Calcula a media das luminancias da imagem
    for (i = 0; i<*ysize; i++) {
        for (j = 0; j<*xsize; j++) {
            acumulador += Image[i][j];
        }
    }
    media = acumulador / ((*xsize)*(*ysize));
    printf("\n Media                         : %d", media);
    //-----------------------------------------


    //Subtrai a média a todos os pixels
    for (i = 0; i<*ysize; i++) {
        for (j = 0; j<*xsize; j++) {
            Image[i][j] -= media;
        }
    }
    //-----------------------------------------


    pointf_out = fopen(outname, "w");
    if (pointf_out == NULL) {
        fprintf(stderr, "Impossivel criar ficheiro de saida: %s\n\n", outname);
        exit(1);
    }

    fprintf(pointf_out, "%d\n", *xsize);
    fprintf(pointf_out, "%d\n", *ysize);
    fprintf(pointf_out, "%d\n", media);

    buffer = 0; /* Bits buffered for ouput		*/
    bits_to_go = 8; /* Number of bits free in buffer	*/
    Bits_Count = 0;

    for (i = 0; i<*ysize; i += G_BsizeY) {
        for (j = 0; j<*xsize; j += G_BsizeX) {
            //Para todos os blocos
            for (i1 = 0; i1 < G_BsizeY; i1++) {
                for (j1 = 0; j1 < G_BsizeX; j1++) {
                    bloco_original[j1 + (i1 * G_BsizeX)] =
                            Image[i + i1][j + j1];
                }
            }

            distorcao = FLT_MAX;
            for (n = 0; n < G_num_codewords; n++) { //Varre todos os elementos do codebook
                aux = quad_err(n);
                if (aux < distorcao) {
                    indice = n;
                    distorcao = aux;
                }
            }

            escreve_indice(indice);

            for (i1 = 0; i1 < G_BsizeY; i1++) {
                for (j1 = 0; j1 < G_BsizeX; j1++) {
                    Image_out[i + i1][j + j1] =
                            G_dic[indice][i1 * G_BsizeX + j1];
                }
            }

        }
    }


    done_outputing_bits(pointf_out);

    end = clock();
    elapsed += ((double) (end - start)) / CLOCKS_PER_SEC;
    start = clock();

    //Soma novamente a média a todos os pixels
    for (i = 0; i<*ysize; i++) {

        for (j = 0; j<*xsize; j++) {
            Image_out[i][j] += media;
        }
    }
    //-----------------------------------------


    psnr = PSNR(Image_orig, Image_out, *ysize, *xsize);
    mse = MSE(Image_orig, Image_out, *ysize, *xsize);
    printf("\n Tempo total de execucao       : %9.3f segundos", elapsed);
    printf("\n Total Bits                    : %ld bits (%ld Bytes)", Bits_Count, Bits_Count / 8);

    taxa = (float) (Bits_Count) / (float) ((*xsize)*(*ysize));

    printf("\n Taxa                          : %.2f bits/pixel", taxa);
    printf("\n PSNR                          : %f dB", psnr);
    printf("\n MSE                           : %f", mse);
    printf("\n-----------------------------------------------------\n\n");


    //write_f_pgm(Image_out, *ysize, *xsize, "Testeout.pgm");

    fclose(pointf_out);
    return EXIT_SUCCESS;
}
//Fim da funcao main

void escreve_indice(int indice) {
    int k;
    int mask;
    int bit_to_write;

    mask = pow(2, G_bits_indice - 1);
    for (k = 0; k < G_bits_indice; k++) {
        bit_to_write = (mask & indice) / mask;
        indice = indice << 1;
        output_bit(bit_to_write, pointf_out);
    }

}




//******************************************************************************
//*                                                                            *
//*     Calcula o erro quadrático entre um vectro do set de treino             *
//*	e um vector do codebook                                                *
//*                                                                            *
//******************************************************************************

void carrega_dicionario(void) {
    int i, j;

    pointf_dic = fopen(dic_name, "r");
    if (pointf_dic == NULL) {
        fprintf(stderr, "Impossivel abrir dicionario: %s\n\n", dic_name);
        exit(1);
    }

    fscanf(pointf_dic, "%d\n", &G_num_codewords);
    fscanf(pointf_dic, "%d\n", &G_BsizeX);
    fscanf(pointf_dic, "%d\n", &G_BsizeY);

    printf("\n-----------------------------------------------------");
    printf("\n Carregou dicionario %s", dic_name);
    printf("\n %d blocos de %dx%d pixels",
            G_num_codewords, G_BsizeY, G_BsizeX);
    printf("\n-----------------------------------------------------");
    fflush(stdout);

    G_dic = intmatrix(G_num_codewords, G_BsizeY * G_BsizeX);

    for (i = 0; i < G_num_codewords; i++) {
        for (j = 0; j < G_BsizeX * G_BsizeY; j++) {
            fscanf(pointf_dic, "%d\t", &G_dic[i][j]);
        }
        (void) fscanf(pointf_dic, "\n");
    }

    fclose(pointf_dic);
}



//******************************************************************************
//*                                                                            *
//*     Calcula o erro quadrático entre um vectro do set de treino             *
//*	e um vector do codebook                                                *
//*                                                                            *
//******************************************************************************

float quad_err(int indice_dic) {
    int i;
    float tmp = 0;

    for (i = 0; i < G_BsizeX * G_BsizeY; i++) {
        tmp += ((G_dic[indice_dic][i] - bloco_original[i]) *
                (G_dic[indice_dic][i] - bloco_original[i]));
    }
    return tmp;
} /* Enf of quad_err */






//******************************************************************************
//*                                                                            *
//*     Apresenta informacao e sintaxe do programa			       *
//*                                                                            *
//******************************************************************************

void ajuda(char *prgname) {
    printf("\n---------------------------------------------------------------------------------\n");
    printf(" Programa de codificação de imagens baseado em quantização vectorial\n");
    printf("---------------------------------------------------------------------------------\n");
    printf(" Sintaxe: %s <ficheiro origem> <ficheiro dicionario> <ficheiro destino>\n", prgname);
    printf("---------------------------------------------------------------------------------\n\n");
}


/************************************************************************************/
/*                                                                                  */
/* READS THE INFORMATION OF A PGM FILE HEADER                                       */
/* usage: read_header_pgm (int *ximg,int *yimg,char *originalfilename);             */
/*                                                                                  */

/************************************************************************************/
void read_header_pgm(int *ysize, int *xsize, char *Fname) {
    FILE *pointf;
#ifdef UNUSED
    char pelimg[1][20];
#endif /*UNUSED */
    char dummy[12];
    char aux;
    int lido = 0;

    pointf = fopen(Fname, "r");
    if (pointf == NULL) {
        fprintf(stderr, "\nCannot open input file: %s\n", Fname);
        exit(1);
    }

    fscanf(pointf, "%s\n", dummy);


    //*******************************
    do {
        aux = getc(pointf);
        //fflush(stdout);
        if (aux == '#') {
            do {
                aux = getc(pointf);
            } while ((aux != '\n') && (aux != EOF));


        } else {
            ungetc(aux, pointf);
            lido = 1;
        }
    } while (lido == 0);

    fscanf(pointf, " %d %d", xsize, ysize);
    //********************************

    fclose(pointf); /* closes file */
}

/* End of read_header_pgm function */




#ifdef CODE_NOT_COMPILED
/*                                                                                  */
/* read_f_pgm: READS A PGM FILE                                                     */
/* usage: read_f_pgm (unsigned char **pelimg, int *ysize, int *xsize, char *Fname)  */
/*                                                                                  */
/* Inputs:                                                                          */
/* Fname - File name                                                                */
/*                                                                                  */
/* Returns 
 * - the image in char **pellimg
   - number of rows 
   - number of columns                                                              */
/*                                                                                  */
#endif /* CODE_NOT_COMPILED */

void read_f_pgm(int **pelimg, int *ysize, int *xsize, char *Fname) {
    int i, j;
    FILE *pointf;
    char dummy[15];
    int aux[1];
    char aux1;
    int lido = 0;


    pointf = fopen(Fname, "r");
    if (pointf == NULL) {
        fprintf(stderr, "\nCannot open input file: %s\n", Fname);
        exit(1);
    }

    fscanf(pointf, "%s\n", dummy);

    //*******************************
    do {
        aux1 = getc(pointf);
        //fflush(stdout);
        if (aux1 == '#') {
            do {
                aux1 = getc(pointf);
            } while ((aux1 != '\n') && (aux1 != EOF));
        } else {
            ungetc(aux1, pointf);
            lido = 1;
        }
    } while (lido == 0);



    fscanf(pointf, " %d %d", xsize, ysize);
    //********************************

    fscanf(pointf, "%d\n", aux);

    for (i = 0; i<*ysize; i++) {
        for (j = 0; j<*xsize; j++) {
            pelimg[i][j] = (unsigned char) fgetc(pointf);
        }
    }

    fclose(pointf); /* closes file */
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
int **intmatrix(int nr, int nc) {
    int i;
    int **m;

    m = (int **) malloc((unsigned) (nr) * sizeof (int *));
    if (!m) {
        printf("intmatrix() - allocation failure 1 \n");
        exit(1);
    }

    for (i = 0; i < nr; i++) {
        m[i] = (int *) malloc((unsigned) (nc) * sizeof (int));
        if (!m[i]) {
            printf("intmatrix() - allocation failure 2 \n");
            exit(1);
        }
    }

    return m;
}

int **intmatrix2(int nr, int nc) {
    int i;
    int vector = (int **) malloc(nr * sizeof (int *));
    for (i = 0; i < nr; i++) {
        // vector[i]=&(bloco[i*nc]);
    }

    // aceder a um bloco
    // bloco [r*nc+c]
    return 0;
}





/************************************************************************************/
/* Peak Signal Noise Ratio                                                          */

/************************************************************************************/
double PSNR(int **origblk, int **cmpblk, int nline, int npixel) {
    int i, j;
    double psnr;

    psnr = 0.;
    for (j = 0; j < nline; j++)
        for (i = 0; i < npixel; i++)
            psnr += ((double) *(origblk[j] + i)-(double) *(cmpblk[j] + i)) * ((double) *(origblk[j] + i)-(double) *(cmpblk[j] + i));

    if (psnr == 0.0) psnr = 0.000000000000000000000000000000000000000000000000000001;
    psnr = 10.0 * log10((double) (RANGEY * RANGEY) / (psnr / (double) (nline * npixel)));

    return psnr;
}
/* End of Psnr function */

/************************************************************************************/
/* Mean Squared Error                                                               */

/************************************************************************************/
double MSE(int **origblk, int **cmpblk, int nline, int npixel) {
    int i, j;
    long cnt = 0;
    double mse;

    mse = 0.;
    for (j = 0; j < nline; j++)
        for (i = 0; i < npixel; i++) {
            mse += ((double) *(origblk[j] + i)-(double) *(cmpblk[j] + i)) * ((double) *(origblk[j] + i)-(double) *(cmpblk[j] + i));
            cnt++;
        }

    return (mse / cnt);
}
/* End of MSE function */

/* INITIALISE BIT OUTPUT */
void start_outputing_bits() {
    buffer = 0; /* Bits buffered for ouput		*/
    bits_to_go = 8; /* Number of bits free in buffer	*/
    Bits_Count = 0;
}

/* OUTPUT A BIT  */
void output_bit(bit, output_file)
int bit;
FILE* output_file;
{
    buffer >>= 1; /* Put bit in top of buffer	*/
    if (bit) buffer |= 0x80;
    bits_to_go -= 1;
    Bits_Count += 1;
    if (bits_to_go == 0) /* Output buffer if it is	*/ { /* now full					*/
        putc(buffer, output_file);
        /*       printf("  %x\n", buffer); */
        bits_to_go = 8;
    }
}

/* FLUSH OUT THE LAST BITS */
void done_outputing_bits(FILE* output_file) {
    putc(buffer >> bits_to_go, output_file);
    //fprintf(stderr," Total Bits: %d bits\n",Bits_Count);
}



/************************************************************************************/
/* WRITE OUTPUT LUMINANCE FILE - PGM format*/

/************************************************************************************/
void write_f_pgm(int **im_matrix, int nline, int npixel, char *filename) {
    int i;
    int pointfo;
    char header_pgm[20];
    int npixel_orig, nline_orig;

    npixel_orig = npixel;
    nline_orig = nline;

    if ((pointfo = creat(filename, PERMS)) == -1) /* creates file */ {
        fprintf(stderr, "Write_f_PGM - cannot create PGM file -> %s", filename);
        exit(1);
    }

    /* Defines the header of file - type pgm */
    *(header_pgm + 0) = 'P';
    *(header_pgm + 1) = '5';
    *(header_pgm + 2) = '\n';

    i = 3;
    if (npixel >= 1000) {
        *(header_pgm + i) = (unsigned char) (npixel / 1000 + 48);
        npixel = npixel % 1000;
        i++;
    }
    *(header_pgm + i) = (unsigned char) (npixel / 100 + 48);
    i++;
    *(header_pgm + i) = (unsigned char) (npixel % 100) / 10 + 48;
    i++;
    *(header_pgm + i) = (unsigned char) ((npixel % 100) % 10 + 48);
    i++;

    *(header_pgm + i) = ' ';
    i++;

    if (nline >= 1000) {
        *(header_pgm + i) = (unsigned char) (nline / 1000 + 48);
        nline = nline % 1000;
        i++;
    }
    *(header_pgm + i) = (unsigned char) (nline / 100 + 48);
    i++;
    *(header_pgm + i) = (unsigned char) ((nline % 100) / 10 + 48);
    i++;
    *(header_pgm + i) = (unsigned char) ((nline % 100) % 10 + 48);
    i++;

    *(header_pgm + i) = '\n';
    i++;

    *(header_pgm + i) = (unsigned char) (RANGE_LUMINANCE / 100 + 48);
    i++;
    *(header_pgm + i) = (unsigned char) ((RANGE_LUMINANCE % 100) / 10 + 48);
    i++;
    *(header_pgm + i) = (unsigned char) ((RANGE_LUMINANCE % 100) % 10 + 48);
    i++;

    *(header_pgm + i) = '\n';
    i++;

    write(pointfo, (char *) header_pgm, i);

    unsigned char **Image_tmp;
    int j;
    Image_tmp = ucmatrix(0, nline, 0, npixel);

    for (i = 0; i < nline; i++)
        for (j = 0; j < npixel; j++)
            Image_tmp[i][j] = (unsigned char) im_matrix[i][j];

    for (i = 0; i < (nline_orig); i++)
        write(pointfo, Image_tmp[i], npixel_orig);

    close(pointfo); /* closes file */
}

/*END of write_f_pgm */

/*************************************************************************************/
/*                                                                                   */
/*  UCMATRIX - Allocates memory for a matrix of variables of type unsigned char      */
/*                                                                                   */
/*  Inputs:                                                                          */
/*    1st and last row / 1st and last column                                         */
/*  Returns a poiter to a unsigned char matrix (unsigned char **)                    */
/*                                                                                   */

/*************************************************************************************/
unsigned char **ucmatrix(int nrl, int nrh, int ncl, int nch) {
    int i;
    unsigned char **m;

    m = (unsigned char **) malloc((unsigned) (nrh - nrl + 1) * sizeof (char*));
    if (!m) {
        printf("ucmatrix() - allocation failure 1 \n");
        exit(1);
    }
    m -= nrl;

    for (i = nrl; i <= nrh; i++) {
        m[i] = (unsigned char *) malloc((unsigned) (nch - ncl + 1) * sizeof (char));
        if (!m[i]) {
            printf("ucmatrix() - allocation failure 2 \n");
            exit(1);
        }
        m[i] -= ncl;
    }
    return m;
}









