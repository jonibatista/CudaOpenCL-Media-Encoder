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


///**************************************************************************///
///*                        CONSTANTS ERROR                                 *///
///**************************************************************************///
#define ERROR_INVALID_PARAMETERS 1
#define ERROR_ALLOCATE_MEMORY 2
#define ERROR_OPEN_FILE
///******************* END CONSTANTES ERROR DEFINITION **********************///


///**************************************************************************///
///*                             CONSTANTS                                  *///
///**************************************************************************///
#define RANGE_LUMINANCE 255 /* Range of image luminance values */
#define PERMS       	0644     /* File acess permits:RW for the users and R for the others */
#define RANGEY      	 255     /* Range level of luminance */

#define Clip1(a)            ((a)>255?255:((a)<0?0:(a)))
///************************ END CONSTANTES DEFINITION ***********************///


///**************************************************************************///
///*                          PROTOTYPES DEFINITION                         *///
///**************************************************************************///
void help(char *prgname);
void read_header_pgm(int *ysize, int *xsize, char *file_name);
void read_f_pgm(int **pelimg, int *ysize, int *xsize, char *file_name);
int **int_matrix(int nr, int nc);
float **floatmatrix(int nr, int nc);
float quad_err(int index_dic, int *block_size_x, int *block_size_y, int *original_block);
void load_dictionary(char *file_name, int *num_codewords, int *block_size_x, int *block_size_y);
double calculate_psnr(int **origblk, int **cmpblk, int nline, int npixel);
double calculate_mse(int **origblk, int **cmpblk, int nline, int npixel);
void write_index(int index, int bits_index, long *bits_count, int *bits_to_go, int *buffer, FILE *pointf_out);

void write_f_pgm(int **im_matrix, int nline, int npixel, char *filename);
unsigned char **ucmatrix(int nrl, int nrh, int ncl, int nch);

//Funcoes de escrita de bits para o ficheiro
void output_bit(int bit, FILE* output_file, int *buffer, int *bits_to_go, long *bits_count);
void done_outputing_bits(FILE* output_file, int *buffer, int *bits_to_go);
///************************ END PROTOTYPES DEFINITION ***********************///


///**************************************************************************///
///*                              GLOBAL VARIABLES                          *///
///**************************************************************************///
int **G_dic;

//******************************************************************************
//*                                                                            *
//*     Funcao main                                                            *
//*                                                                            *
//*                                                                            *
//******************************************************************************

int main(int argc, char *argv[]) {
    int **image, **image_orig, **image_out;
    int *original_block;

    int i, j, i1, j1, n;
    int index = -1; /* Dummy value -- index is set before its usage */
    int average;
    float distortion;
    float aux = 0;
    double psnr, mse;

    clock_t start, end;
    double elapsed = 0;
    start = clock();

    /* THE BIT BUFFER */
    int buffer = 0; /* Bits buffered for output 		*/
    int bits_to_go = 8; /* Number of bits still in buffer	        */
    long bits_count = 0;

    /*
     * The following variables are loaded from the dictionary  file
     **/
    int block_size_x = 0;
    int block_size_y = 0;
    int num_codewords = 0;

    int bits_index = 0;

    int ysize = 0, xsize = 0; /* The dimensions of the original image */

    char *inname, *outname, *dic_name; //Nome dos ficheiros de input e output

    FILE *pointf_out;

    //Verifica argumentos
    if (argc < 3) {
        help(argv[0]);
        return 0;
    }

    inname = argv[1];
    dic_name = argv[2];
    outname = argv[3];


    //Carrega dicionario
    load_dictionary(dic_name, &num_codewords, &block_size_x, &block_size_y);
    bits_index = ceil(log(num_codewords) / log(2));

    original_block = (int *) calloc(block_size_x*block_size_y, sizeof (int));
    if (!original_block) {
        printf("int_matrix() - allocation failure 1 \n");
        exit(1);
    }

    //Le imagem a comprimir
    printf("\n imagem a comprimir            : %s", inname);
    read_header_pgm(&ysize, &xsize, inname); /* Reads the PGM file and returns the picture size */
    image_orig = int_matrix(ysize, xsize);
    image = int_matrix(ysize, xsize);
    image_out = int_matrix(ysize, xsize);
    read_f_pgm(image, &ysize, &xsize, inname); /* Reads the PGM file and stores the image in pely */
    printf("\n Tamanho (%dx%d)             : %d pixels", xsize, ysize, xsize * ysize);

    for (i = 0; i < ysize; i++) {
        for (j = 0; j < xsize; j++) {
            image_orig[i][j] = image[i][j];
        }
    }

    int num_blocks = (ysize / block_size_y)*(xsize / block_size_x);
    printf("\n Total de blocos %dx%d na imagem : %d blocos",
            block_size_y, block_size_x, num_blocks);
    printf("\n-----------------------------------------------------");

    //Calcula a average das luminancias da imagem

    for (i = 0; i < ysize; i++) {
        for (j = 0; j < xsize; j++) {
            aux += image[i][j];
        }
    }
    average = aux / ((xsize)*(ysize));
    printf("\n average                         : %d", average);
    fflush(stdout);
    //-----------------------------------------


    //Subtrai a média a todos os pixels
    for (i = 0; i < ysize; i++) {
        for (j = 0; j < xsize; j++) {
            image[i][j] -= average;
        }
    }
    //-----------------------------------------


    pointf_out = fopen(outname, "w");
    if (pointf_out == NULL) {
        fprintf(stderr, "Impossivel criar ficheiro de saida: %s\n\n", outname);
        exit(1);
    }

    fprintf(pointf_out, "%d\n", xsize);
    fprintf(pointf_out, "%d\n", ysize);
    fprintf(pointf_out, "%d\n", average);

    for (i = 0; i < ysize; i += block_size_y) {
        for (j = 0; j < xsize; j += block_size_x) {
            //Para todos os blocos
            for (i1 = 0; i1 < block_size_y; i1++) {
                for (j1 = 0; j1 < block_size_x; j1++) {
                    original_block[j1 + (i1 * block_size_x)] =
                            image[i + i1][j + j1];
                }
            }

            distortion = FLT_MAX;
            for (n = 0; n < num_codewords; n++) { //Varre todos os elementos do codebook
                aux = quad_err(n, &block_size_x, &block_size_y, original_block);
                if (aux < distortion) {
                    index = n;
                    distortion = aux;
                }
            }

            write_index(index, bits_index, &bits_count, &bits_to_go, &buffer, pointf_out);

            for (i1 = 0; i1 < block_size_y; i1++) {
                for (j1 = 0; j1 < block_size_x; j1++) {
                    image_out[i + i1][j + j1] =
                            G_dic[index][i1 * block_size_x + j1];
                }
            }

        }
    }


    done_outputing_bits(pointf_out, &buffer, &bits_to_go);

    end = clock();
    elapsed += ((double) (end - start)) / CLOCKS_PER_SEC;
    start = clock();

    //Soma novamente a média a todos os pixels
    for (i = 0; i < ysize; i++) {

        for (j = 0; j < xsize; j++) {
            image_out[i][j] += average;
        }
    }
    //-----------------------------------------


    psnr = calculate_psnr(image_orig, image_out, ysize, xsize);
    mse = calculate_mse(image_orig, image_out, ysize, xsize);
    printf("\n Tempo total de execucao       : %9.3f segundos", elapsed);
    printf("\n Total Bits                    : %ld bits (%ld Bytes)", bits_count, bits_count / 8);

    float rate = (float) (bits_count) / (float) (xsize * ysize);

    printf("\n Taxa                          : %.2f bits/pixel", rate);
    printf("\n psnr                          : %f dB", psnr);
    printf("\n mse                           : %f", mse);
    printf("\n-----------------------------------------------------\n\n");


    //write_f_pgm(image_out, *ysize, *xsize, "Testeout.pgm");

    fclose(pointf_out);
    return EXIT_SUCCESS;
}
//Fim da funcao main

void write_index(int index, int bits_index, long *bits_count, int *bits_to_go, int *buffer, FILE *pointf_out) {
    int k;
    int mask;
    int bit_to_write;

    mask = pow(2, bits_index - 1);
    for (k = 0; k < bits_index; k++) {
        bit_to_write = (mask & index) / mask;
        index = index << 1;
        output_bit(bit_to_write, pointf_out, buffer, bits_to_go, bits_count);
    }

}




//******************************************************************************
//*                                                                            *
//*     Calcula o erro quadrático entre um vectro do set de treino             *
//*	e um vector do codebook                                                *
//*                                                                            *
//******************************************************************************

void load_dictionary(char *file_name, int *num_codewords, int *block_size_x, int *block_size_y) {
    int i, j;
    FILE *pointf_dic;

    pointf_dic = fopen(file_name, "r");
    if (pointf_dic == NULL) {
        fprintf(stderr, "Impossivel abrir dicionario: %s\n\n", file_name);
        exit(1);
    }

    fscanf(pointf_dic, "%d\n", num_codewords);
    fscanf(pointf_dic, "%d\n", block_size_x);
    fscanf(pointf_dic, "%d\n", block_size_y);

    printf("\n-----------------------------------------------------");
    printf("\n Carregou dicionario %s", file_name);
    printf("\n %d blocos de %dx%d pixels",
            *num_codewords, *block_size_y, *block_size_x);
    printf("\n-----------------------------------------------------");
    fflush(stdout);

    G_dic = int_matrix(*num_codewords, *block_size_y * (*block_size_x));

    for (i = 0; i<*num_codewords; i++) {
        for (j = 0; j < *block_size_x * (*block_size_y); j++) {
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

float quad_err(int index_dic, int *block_size_x, int *block_size_y, int *original_block) {
    int i;
    float tmp = 0;

    for (i = 0; i < *block_size_x * (*block_size_y); i++) {
        tmp += ((G_dic[index_dic][i] - original_block[i]) *
                (G_dic[index_dic][i] - original_block[i]));
    }
    return tmp;
} /* Enf of quad_err */






//******************************************************************************
//*                                                                            *
//*     Apresenta informacao e sintaxe do programa			       *
//*                                                                            *
//******************************************************************************

void help(char *prgname) {
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
/**
 *  <p> Reads the information of a pgm file to calculate the horizontal and vertical size.</p>
 * 
 * @param ysize image vertical dimension
 * @param xsize image horizontal dimensio
 * @param file_name file name of the image that will be coded
 */
void read_header_pgm(int *ysize, int *xsize, char *file_name) {
    FILE *pointf;
    char dummy[12];
    char aux;
    int lido = 0;

    pointf = fopen(file_name, "r");
    if (pointf == NULL) {
        fprintf(stderr, "\nCannot open input file: %s\n", file_name);
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


/**
 *
 * <p> Copy the image to memory </p>
 *
 * @param pelimg vector where the image will be saved
 * @param ysize image vertical dimension
 * @param xsize image horizontal dimensio
 * @param file_name file name of the image that will be coded
 */
void read_f_pgm(int **pelimg, int *ysize, int *xsize, char *file_name) {
    int i, j;
    FILE *pointf;
    char dummy[15];
    int aux[1];
    char aux1;
    int lido = 0;


    pointf = fopen(file_name, "r");
    if (pointf == NULL) {
        fprintf(stderr, "\nCannot open input file: %s\n", file_name);
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
/*  int_matrix - Allocates memory for a matrix of variables of type int               */
/*                                                                                   */
/*  Inputs:                                                                          */
/*    number of rows / nunmber of columnc                                            */
/*  Returns a poiter to a int matrix (int **)                                        */
/*                                                                                   */
/*************************************************************************************/
int **int_matrix(int nr, int nc) {
    int i;
    int **m;

    m = (int **) malloc((unsigned) (nr) * sizeof (int *));
    if (!m) {
        printf("int_matrix() - allocation failure 1 \n");
        exit(1);
    }

    for (i = 0; i < nr; i++) {
        m[i] = (int *) malloc((unsigned) (nc) * sizeof (int));
        if (!m[i]) {
            printf("int_matrix() - allocation failure 2 \n");
            exit(1);
        }
    }

    return m;
}

float **floatmatrix(int nr, int nc) {
    int i;
    float **m;

    m = (float **) malloc(nr * nc * sizeof (float *));
    if (!m) {
        printf("floatmatrix() - allocation failure 1 \n");
        exit(1);
    }

    for (i = 0; i < nr; i++) {
        m[i] = (float *) calloc(nc, sizeof (float));
        if (!m[i]) {
            printf("floatmatrix() - allocation failure 2 \n");
            exit(1);
        }
    }

    return m;
}


/************************************************************************************/
/* Peak Signal Noise Ratio                                                          */

/************************************************************************************/
double calculate_psnr(int **origblk, int **cmpblk, int nline, int npixel) {
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
/* End of psnr function */

/************************************************************************************/
/* Mean Squared Error                                                               */

/************************************************************************************/
double calculate_mse(int **origblk, int **cmpblk, int nline, int npixel) {
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
/* End of mse function */

/* OUTPUT A BIT  */
void output_bit(int bit, FILE* output_file, int *buffer, int *bits_to_go, long *bits_count) {
    *buffer >>= 1; /* Put bit in top of buffer	*/
    if (bit) *buffer |= 0x80;
    *bits_to_go -= 1;
    *bits_count += 1;
    if (*bits_to_go == 0) /* Output buffer if it is	*/ { /* now full					*/
        putc(*buffer, output_file);
        /*       printf("  %x\n", buffer); */
        *bits_to_go = 8;
    }
}

/* FLUSH OUT THE LAST BITS */
void done_outputing_bits(FILE* output_file, int *buffer, int *bits_to_go) {
    putc(*buffer >> *bits_to_go, output_file);
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

    unsigned char **image_tmp;
    int j;
    image_tmp = ucmatrix(0, nline, 0, npixel);

    for (i = 0; i < nline; i++)
        for (j = 0; j < npixel; j++)
            image_tmp[i][j] = (unsigned char) im_matrix[i][j];

    for (i = 0; i < (nline_orig); i++)
        write(pointfo, image_tmp[i], npixel_orig);

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





