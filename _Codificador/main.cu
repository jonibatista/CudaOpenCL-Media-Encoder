////////////////////////////////////////////////////////////////////////////////
///   Copyright (C) 2008 by Nelson Carreira Francisco                        ///
///   eng.nelsito@gmail.com                                                  ///
///                                                                          ///
///   This program is free software; you can redistribute it and/or modify   ///
///   it under the terms of the GNU General Public License as published by   ///
///   the Free Software Foundation; either version 2 of the License, or      ///
///   (at your option) any later version.                                    ///
///                                                                          ///
///   This program is distributed in the hope that it will be useful,        ///
///   but WITHOUT ANY WARRANTY; without even the implied warranty of         ///
///   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ///
///   GNU General Public License for more details.                           ///
///                                                                          ///
///   You should have received a copy of the GNU General Public License      ///
///   along with this program; if not, write to the                          ///
///   Free Software Foundation, Inc.,                                        ///
///   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              ///
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
///   Implementacaoo de um codificador de imagens baseado em                 ///
///   Quantificacao vectorial                                                ///
///   Nelson Carreira Francisco                                              ///
///   Optimized by César Ferreira e Jóni Batista                           ///
////////////////////////////////////////////////////////////////////////////////


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
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "params.h"


////////////////////////////////////////////////////////////////////////////////
///                         CONSTANTS ERROR                                  ///
////////////////////////////////////////////////////////////////////////////////
#define ERROR_INVALID_PARAMETERS 1
#define ERROR_ALLOCATE_MEMORY 2
#define ERROR_OPEN_FILE
#define ERROR_INVALID_PARAMETERS 1



////////////////////////////////////////////////////////////////////////////////
///                              CONSTANTS                                   ///
////////////////////////////////////////////////////////////////////////////////
#define RANGE_LUMINANCE 255	// Range of image luminance values
#define PERMS       	0644	// File acess permits:RW for the users and R for the others
#define RANGEY      	 255	// Range level of luminance

#define Clip1(a)            ((a)>255?255:((a)<0?0:(a)))
/**
 * <p> Function to process CUDA errors </p>
 *
 * @param err [IN] CUDA error to process (usually the code returned by the cuda function)
 * @param line [IN] line of source code where function is called
 * @param file [IN] name of source file where function is called
 * @return on error, the function terminates the process with EXIT_FAILURE code.
 *
 * source: "CUDA by Example: An Introduction to General-Purpose "
 * GPU Programming", Jason Sanders, Edward Kandrot, NVIDIA, July 2010
 * @note: the function should be called through the macro 'HANDLE_ERROR'
 **/
static void
HandleError (cudaError_t err, const char *file, int line)
{
  if (err != cudaSuccess)
    {
      printf ("[ERROR] '%s' (%d) in '%s' at line '%d'\n",
	      cudaGetErrorString (err), err, file, line);
      exit (EXIT_FAILURE);
    }
}

/**
 * <p>HANDLE_ERROR macro.</p>
 *
 * Wrapping macro for HandleError function (provides "file" and "line" parameters).
 *
 * @param err [IN] CUDA error
 * @return on error, the calling process is terminated
 **/
#define HANDLE_ERROR(err) (HandleError((err), __FILE__, __LINE__ ))


const int G_BlocksPerGrid = 65535;	//65535; => tesla max
const int G_ThreadsPerBlock = 64;	//1024; => tesla max

////////////////////////////////////////////////////////////////////////////////
///                             CUDA KERNEL                                  ///
////////////////////////////////////////////////////////////////////////////////
/**
 * <p>calculate all quad error and create the pgm encoded</p>
 * 
 * @param num_codewords number of rows of the dictionary
 * @param block_size size of the block (horizontal * verticar coordinates)
 * @param dev_dict ditionary data in memory
 * @param dev_pgm pgm image in memory
 * @param dev_pgm_coded encoding result
 */
__global__ void
encoding_pgm (int num_codewords, int pgm_block_size, int *dev_dict,
	      int *dev_pgm, int *dev_pgm_coded, int num_blocks)
{

  if (blockIdx.x < num_blocks)
    {
      __shared__ float cache_err[G_ThreadsPerBlock];
      __shared__ int cache_idx[G_ThreadsPerBlock];

      int i, idx_dict, idx_block;
      int tid = threadIdx.x + (blockIdx.x * blockDim.x);
      long int jump = blockDim.x * gridDim.x;
      int global_size = blockDim.x * gridDim.x - 1;	// lets us know when starts the local stride
      float temp = 0.0;

      if (threadIdx.x == 0)
	{
	  for (i = 0; i < blockDim.x; i++)
	    {
	      cache_err[i] = FLT_MAX;
	      cache_idx[i] = 0;
	    }
	}

      __syncthreads ();

// make the dictionary stride  
      while (tid < (gridDim.x * num_codewords))

	{
	  if (tid > global_size)
	    {
	      idx_dict += threadIdx.x * pgm_block_size;
	    }
	  else
	    {
	      idx_dict = threadIdx.x * pgm_block_size;
	    }
	  i = 0;
	  temp = 0.0;
	  idx_block = 0;

	  while (i < pgm_block_size)
	    {
	      idx_block = (blockIdx.x * pgm_block_size) + i;
	      temp +=
		((dev_dict[idx_dict + i] -
		  dev_pgm[idx_block]) * (dev_dict[idx_dict + i] -
					 dev_pgm[idx_block]));
	      i++;
	    }

	  if (cache_err[threadIdx.x] > temp)
	    {
	      cache_err[threadIdx.x] = temp;
	      cache_idx[threadIdx.x] = idx_dict / pgm_block_size;
	    }
	  tid += jump;
	}

      __syncthreads ();

      if (threadIdx.x == 0)
	{
	  float aux = FLT_MAX;

	  for (i = 0; i < blockDim.x; i++)
	    {
	      if (cache_err[i] < aux)
		{
		  aux = cache_err[i];
		  dev_pgm_coded[blockIdx.x] = cache_idx[i];
		}
	    }
	}
    }
}


////////////////////////////////////////////////////////////////////////////////
///                           PROTOTYPES DEFINITION                          ///
////////////////////////////////////////////////////////////////////////////////

void read_header_pgm (int *ysize, int *xsize, char *file_name);
void read_file_pgm (int **pelimg, int *ysize, int *xsize, char *file_name);
void v_read_file_pgm (int *pelimg, int *ysize, int *xsize, char *file_name);
int **int_matrix (int nr, int nc);
float **floatmatrix (int nr, int nc);
int *int_vector (int nr, int nc);
void sort_pgm_blocks (int *v_pgm, int *v_sort_pgm, int block_x, int block_y,
		      int xsize, int ysize);
float quad_err (int index_dic, int block_size, int *original_block);
void load_dictionary (char *file_name, int *num_codewords, int *block_size_x,
		      int *block_size_y);
double calculate_psnr (int **origblk, int **cmpblk, int nline, int npixel);
double calculate_mse (int **origblk, int **cmpblk, int nline, int npixel);
void write_index (int index, int bits_index, long *bits_count,
		  int *bits_to_go, int *buffer, FILE * pointf_out);

void write_f_pgm (int **im_matrix, int nline, int npixel, char *filename);
unsigned char **ucmatrix (int nrl, int nrh, int ncl, int nch);

//Funcoes de escrita de bits para o ficheiro
void output_bit (int bit, FILE * output_file, int *buffer, int *bits_to_go,
		 long *bits_count);
void done_outputing_bits (FILE * output_file, int *buffer, int *bits_to_go);
int max_number_threads ();
int max_number_blocks ();

// CUDA calls
void kernel_threads_per_codwords (int block_size, int num_codewords,
				  int xsize, int ysize, int *v_pgm_sorted,
				  int *v_pgm_coded);


////////////////////////////////////////////////////////////////////////////////
///                               GLOBAL VARIABLES                           ///
////////////////////////////////////////////////////////////////////////////////
int *G_dic;


////////////////////////////////////////////////////////////////////////////////
///                                  FUNCTIONS                               ///
////////////////////////////////////////////////////////////////////////////////

int
main (int argc, char *argv[])
{
  int **image_orig, **image_out;
  int *v_pgm, *v_pgm_sorted, *v_pgm_coded;

  struct gengetopt_args_info args_info;

  int *original_block;

  int i, j, i1, j1;
  int index = -1;		/* Dummy value -- index is set before its usage */
  int average;
  float aux = 0;
  double psnr, mse;

  clock_t start, end;
  double elapsed = 0;
  start = clock ();

  /* THE BIT BUFFER */
  int buffer = 0;		/* Bits buffered for output                 */
  int bits_to_go = 8;		/* Number of bits still in buffer               */
  long bits_count = 0;

  /*
   * The following variables are loaded from the dictionary  file
   **/
  int block_size = 0;
  int block_size_x = 0;
  int block_size_y = 0;
  int num_codewords = 0;

  int bits_index = 0;

  int ysize = 0, xsize = 0;	/* The dimensions of the original image */

  char *inname, *outname, *dic_name;	//Nome dos ficheiros de input e output

  FILE *pointf_out;


  printf ("=====================================================");
  printf ("\n Total de blocos do GPU: %d", G_BlocksPerGrid);
  printf ("\n Total de Threads de cada Bloco do GPU:  %d", G_ThreadsPerBlock);
  printf ("\n=====================================================");

  // validate parameters
  if (cmdline_parser (argc, argv, &args_info) != 0)
    {
      exit (ERROR_INVALID_PARAMETERS);
    }

  inname = args_info.image_arg;
  dic_name = args_info.dictionary_arg;
  outname = args_info.file_arg;

  //Carrega dicionario
  load_dictionary (dic_name, &num_codewords, &block_size_x, &block_size_y);
  bits_index = ceil (log (num_codewords) / log (2));

  block_size = block_size_x * block_size_y;
  int dict_statistics[num_codewords][block_size];
  for (i = 0; i < num_codewords; i++)
    for (j = 0; j < block_size; j++)
      dict_statistics[i][j] = G_dic[i * block_size + j];

  original_block = (int *) calloc (block_size, sizeof (int));
  if (!original_block)
    {
      printf ("int_matrix() - allocation failure 1 \n");
      exit (1);
    }
 
  //Le imagem a comprimir
  printf ("\n imagem a comprimir            : %s", inname);
  read_header_pgm (&ysize, &xsize, inname);	/* Reads the PGM file and returns the picture size */

  // load pgm to vector
  v_pgm = int_vector (ysize, xsize);
  v_read_file_pgm (v_pgm, &ysize, &xsize, inname);

  // create the vector that will contain the coded pgm
  v_pgm_coded = int_vector (ysize / block_size_y, xsize / block_size_x);

  // old stuff (using matrixes)
  image_orig = int_matrix (ysize, xsize);
  image_out = int_matrix (ysize, xsize);


  printf ("\n Tamanho (%dx%d)             : %d pixels", xsize, ysize,
	  xsize * ysize);

  for (i = 0; i < ysize; i++)
    {
      for (j = 0; j < xsize; j++)
	{
	  image_orig[i][j] = v_pgm[i * xsize + j];
	}
    }

  int num_blocks = (ysize / block_size_y) * (xsize / block_size_x);
  printf ("\n Total de blocos %dx%d na imagem : %d blocos",
	  block_size_y, block_size_x, num_blocks);
  printf ("\n-----------------------------------------------------");

  //Calcula a average das luminancias da imagem

  for (i = 0; i < ysize; i++)
    {
      for (j = 0; j < xsize; j++)
	{
	  aux += v_pgm[i * xsize + j];
	}
    }
  average = aux / ((xsize) * (ysize));
  printf ("\n average                         : %d", average);
  fflush (stdout);
  //-----------------------------------------


  //Subtrai a média a todos os pixels
  for (i = 0; i < ysize; i++)
    {
      for (j = 0; j < xsize; j++)
	{
	  v_pgm[i * xsize + j] -= average;
	}
    }
  //-----------------------------------------


  pointf_out = fopen (outname, "w");
  if (pointf_out == NULL)
    {
      fprintf (stderr, "Impossivel criar ficheiro de saida: %s\n\n", outname);
      exit (1);
    }


  fprintf (pointf_out, "%d\n", xsize);
  fprintf (pointf_out, "%d\n", ysize);
  fprintf (pointf_out, "%d\n", average);

  // sort the pixels of the each block
  v_pgm_sorted = int_vector (ysize, xsize);
  sort_pgm_blocks (v_pgm, v_pgm_sorted, block_size_x, block_size_y, xsize,
		   ysize);

  //
  // CUDA STUFF
  //

  //call encoding function 
  kernel_threads_per_codwords (block_size, num_codewords, xsize, ysize,
			       v_pgm_sorted, v_pgm_coded);

//
//END CUDA STUFF
//

// compare gpu encoded with cpu encoded
/*  for (i = 0; i < xsize * ysize / block_size; i++)
    {
      printf ("%d ", v_pgm_coded[i]);
    }
*/

  // estatistics
  for (i = 0; i < ysize; i += block_size_y)
    {
      for (j = 0; j < xsize; j += block_size_x)
	{
	  index =
	    v_pgm_coded[(i / block_size_y) * (xsize / block_size_x) +
			(j / block_size_x)];
	  for (i1 = 0; i1 < block_size_y; i1++)
	    {
	      for (j1 = 0; j1 < block_size_x; j1++)
		{
		  image_out[i + i1][j + j1] =
		    dict_statistics[index][i1 * block_size_x + j1];
		}
	    }

	}
    }

  // write the coded pgm to file
  for (i = 0; i < (ysize / block_size_y); i++)
    {
      for (j = 0; j < (xsize / block_size_x); j++)
	{
	  write_index (v_pgm_coded[i * (xsize / block_size_x) + j],
		       bits_index, &bits_count, &bits_to_go, &buffer,
		       pointf_out);
	}
    }


  done_outputing_bits (pointf_out, &buffer, &bits_to_go);

  end = clock ();
  elapsed += ((double) (end - start)) / CLOCKS_PER_SEC;
  start = clock ();

  //Soma novamente a média a todos os pixels
  for (i = 0; i < ysize; i++)
    {

      for (j = 0; j < xsize; j++)
	{
	  image_out[i][j] += average;
	}
    }
  //-----------------------------------------


  psnr = calculate_psnr (image_orig, image_out, ysize, xsize);
  mse = calculate_mse (image_orig, image_out, ysize, xsize);
  printf ("\n Tempo total de execucao       : %9.3f segundos", elapsed);
  printf ("\n Total Bits                    : %ld bits (%ld Bytes)",
	  bits_count, bits_count / 8);

  float rate = (float) (bits_count) / (float) (xsize * ysize);

  printf ("\n Taxa                          : %.2f bits/pixel", rate);
  printf ("\n psnr                          : %f dB", psnr);
  printf ("\n mse                           : %f", mse);
  printf ("\n-----------------------------------------------------\n\n");


  fclose (pointf_out);

  // now free the memory
  for (i = 0; i < ysize; i++)
    {
      free (image_orig[i]);
      free (image_out[i]);
    }

  //free memory
  free (image_orig);
  image_orig = NULL;
  free (image_out);
  image_out = NULL;
  free (v_pgm_coded);
  v_pgm_coded = NULL;
  free (v_pgm);
  v_pgm = NULL;
  free (v_pgm_sorted);
  v_pgm_sorted = NULL;
  free (G_dic);
  G_dic = NULL;

  return EXIT_SUCCESS;
}



/**
 * <p> </p>
 *
 * @param block_size the image block size
 * @param num_codewords the number of rows in the dictionary file
 * @param xsize horizontal image size
 * @param ysize vertical image size
 * @param v_pgm_sorted the vector with the image src pixeis
 * @param v_pgm_coded the vector with the encoded image (dictionary index)
 */
void
kernel_threads_per_codwords (int block_size, int num_codewords, int xsize,
			     int ysize, int *v_pgm_sorted, int *v_pgm_coded)
{
  int i, i1, j1;
  // the gpu device vectors
  int *dev_pgm;
  int *dev_dict;
  int *dev_pgm_coded;


// the temporary vector of one parte of the image.
// because the number of GPU blocks could be lesse than the number of image blocks
  int *v_pgm_temp, *v_pgm_coded_temp;

// iterator for gpu calls. 
  int it_gpu_results;

  // alloc memory for temporary vectors
  v_pgm_temp =
    (int *) malloc ((unsigned) G_BlocksPerGrid * block_size * sizeof (int *));
  if (!v_pgm_temp)
    {
      printf ("int_vector() - allocation failure 1 \n");
      exit (1);
    }

  v_pgm_coded_temp =
    (int *) malloc ((unsigned) G_BlocksPerGrid * sizeof (int *));
  if (!v_pgm_coded_temp)
    {
      printf ("int_vector() - allocation failure 1 \n");
      exit (1);
    }
  for (i = 0; i < G_BlocksPerGrid; i++)
    v_pgm_coded_temp[i] = 0;
  int num_gpu_calls = ((xsize * ysize) / block_size) / G_BlocksPerGrid;
  int num_blocks_grid = G_BlocksPerGrid, missing_blocks = 0;

  if (((xsize * ysize) / block_size) % G_BlocksPerGrid != 0)
    {
      missing_blocks =
	((xsize * ysize) -
	 (G_BlocksPerGrid * block_size * num_gpu_calls)) / block_size;
    }


  // alloc memory to cuda vectors
  HANDLE_ERROR (cudaMalloc
		((void **) &dev_pgm,
		 (size_t) G_BlocksPerGrid * block_size * sizeof (int)));
  HANDLE_ERROR (cudaMalloc
		((void **) &dev_dict,
		 (size_t) block_size * num_codewords * sizeof (int)));
  HANDLE_ERROR (cudaMalloc
		((void **) &dev_pgm_coded,
		 (size_t) G_BlocksPerGrid * sizeof (int)));


  HANDLE_ERROR (cudaMemcpy
		(dev_dict, G_dic, num_codewords * block_size * sizeof (int),
		 cudaMemcpyHostToDevice));

  for (i = 0, it_gpu_results = 0;
       i <
       (G_BlocksPerGrid * block_size * num_gpu_calls) +
       (missing_blocks * block_size);
       i += (G_BlocksPerGrid * block_size), it_gpu_results++)
    {

      // copy the images blocks to temp...
      for (i1 = i, j1 = 0; j1 < num_blocks_grid * block_size; i1++, j1++)
	v_pgm_temp[j1] = v_pgm_sorted[i1];

      // copy the vectors data fom host to gpu device
      HANDLE_ERROR (cudaMemcpy
		    (dev_pgm, v_pgm_temp,
		     G_BlocksPerGrid * block_size * sizeof (int),
		     cudaMemcpyHostToDevice));

      // execute GPU KERNEL
      encoding_pgm << <G_BlocksPerGrid, G_ThreadsPerBlock >> >(num_codewords,
							       block_size,
							       dev_dict,
							       dev_pgm,
							       dev_pgm_coded, num_blocks_grid);

      cudaThreadSynchronize ();

      // copy the vector with the pgm coded from dpu decive to host
      HANDLE_ERROR (cudaMemcpy
		    (v_pgm_coded_temp, dev_pgm_coded,
		     (G_BlocksPerGrid) * sizeof (int),
		     cudaMemcpyDeviceToHost));

      //add the short vector result to final vector result
      for (i1 = 0; i1 < num_blocks_grid; i1++)
	{
	  v_pgm_coded[(G_BlocksPerGrid * it_gpu_results) + i1] =
	    v_pgm_coded_temp[i1];
	}

      if (missing_blocks != 0 && num_gpu_calls == it_gpu_results + 1)
	{
	  // there're missing blocks...
	  // this time the gpu will use less blocks
	  num_blocks_grid = missing_blocks;
	}
    }

/*
for(int x=0; x<(xsize * ysize)/block_size; x++){
	printf("\n%d = %d =>", v_pgm_coded2[x], v_pgm_coded[x]);
	if(v_pgm_coded2[x] != v_pgm_coded[x])
               printf("PROBLEM!!!!!!!!!!!!!!!!!!!!");
 }
*/

// cuda free memory
  cudaFree (dev_pgm);
  dev_pgm = NULL;
  cudaFree (dev_dict);
  dev_dict = NULL;
  cudaFree (dev_pgm_coded);
  dev_pgm_coded = NULL;
  free (v_pgm_temp);
  v_pgm_temp = NULL;
  free (v_pgm_coded_temp);
  v_pgm_coded_temp = NULL;

}


/**
 * <p> Gets the max number of threads for this device </p>
 *
 * @return max number of threads for this device
 */
int
max_number_threads ()
{
  cudaDeviceProp prop;
  cudaGetDeviceProperties (&prop, 0);
  return prop.maxThreadsPerBlock;
}


/**
 * <p> Gets the max number of GPU bolcks for this device </p>
 *
 * @return max number of GPU blocks for this device
 */
int
max_number_blocks ()
{
  cudaDeviceProp prop;
  cudaGetDeviceProperties (&prop, 0);
  return prop.maxGridSize[0];
}


/**
 *
 * @param index
 * @param bits_index
 * @param bits_count
 * @param bits_to_go
 * @param buffer
 * @param pointf_out
 */
void
write_index (int index, int bits_index, long *bits_count, int *bits_to_go,
	     int *buffer, FILE * pointf_out)
{
  int k;
  int mask;
  int bit_to_write;

  mask = pow (2, bits_index - 1);
  for (k = 0; k < bits_index; k++)
    {
      bit_to_write = (mask & index) / mask;
      index = index << 1;
      output_bit (bit_to_write, pointf_out, buffer, bits_to_go, bits_count);
    }

}

/**
 * <p> Load the dictionary file to memory. </p>
 *
 * @param file_name the name of the dictionary file
 * @param num_codewords number of blocks of the dictionary
 * @param block_size_x horizontal size of the block
 * @param block_size_y vertical size of the block
 *
 */
void
load_dictionary (char *file_name, int *num_codewords, int *block_size_x,
		 int *block_size_y)
{
  int i, j;
  FILE *pointf_dic;

  pointf_dic = fopen (file_name, "r");
  if (pointf_dic == NULL)
    {
      fprintf (stderr, "Impossivel abrir dicionario: %s\n\n", file_name);
      exit (1);
    }

  if (fscanf (pointf_dic, "%d\n", num_codewords) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }
  if (fscanf (pointf_dic, "%d\n", block_size_x) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }
  if (fscanf (pointf_dic, "%d\n", block_size_y) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }

  printf ("\n-----------------------------------------------------");
  printf ("\n Carregou dicionario %s", file_name);
  printf ("\n %d blocos de %dx%d pixels",
	  *num_codewords, *block_size_y, *block_size_x);
  printf ("\n-----------------------------------------------------");
  fflush (stdout);

  G_dic = int_vector (*num_codewords, *block_size_y * (*block_size_x));

  for (i = 0; i < *num_codewords; i++)
    {
      for (j = 0; j < *block_size_x * (*block_size_y); j++)
	{
	  if (fscanf (pointf_dic, "%d\t",
		      &G_dic[i * (*block_size_x * (*block_size_y)) + j]) ==
	      EOF)
	    {
	      printf ("\nEOF ERROR\n");

	    }

	}
      //(void) fscanf(pointf_dic, "\n");
      if (fscanf (pointf_dic, "\n") == EOF)
	{
	  printf ("\nEOF ERROR\n");

	}
    }

  fclose (pointf_dic);
}

/**
 * <p> Calculate the square error between a vector and a training set of the codebook vector. </p>
 *
 * @param index_dic index of the dictionary row
 * @param block_size size of the block (horizontal coordinate * vertical coordinate)
 * @param original_block the current block
 * @return the square error value
 */
float
quad_err (int index_dic, int block_size, int *original_block)
{
  int i;
  float tmp = 0;

  for (i = 0; i < block_size; i++)
    {
      tmp +=
	((G_dic[index_dic * block_size + i] -
	  original_block[i]) * (G_dic[index_dic * block_size + i] -
				original_block[i]));
    }
  return tmp;
}

/**
 *  <p> Reads the information of a pgm file to calculate the horizontal and vertical size.</p>
 *
 * @param ysize image vertical dimension
 * @param xsize image horizontal dimensio
 * @param file_name file name of the image that will be coded
 */
void
read_header_pgm (int *ysize, int *xsize, char *file_name)
{
  FILE *pointf;
  char dummy[12];
  char aux;
  int read = 0;

  pointf = fopen (file_name, "r");
  if (pointf == NULL)
    {
      fprintf (stderr, "\nCannot open input file: %s\n", file_name);
      exit (1);
    }

  if (fscanf (pointf, "%s\n", dummy) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }


  //*******************************
  do
    {
      aux = getc (pointf);
      //fflush(stdout);
      if (aux == '#')
	{
	  do
	    {
	      aux = getc (pointf);
	    }
	  while ((aux != '\n') && (aux != EOF));


	}
      else
	{
	  ungetc (aux, pointf);
	  read = 1;
	}
    }
  while (read == 0);

  if (fscanf (pointf, " %d %d", xsize, ysize) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }
  //********************************

  fclose (pointf);		/* closes file */
}

/**
 *
 * <p> Copy the image to memory </p>
 *
 * @param pelimg vector where the image will be saved
 * @param ysize image vertical dimension
 * @param xsize image horizontal dimensio
 * @param file_name file name of the image that will be coded
 */
void
read_file_pgm (int **pelimg, int *ysize, int *xsize, char *file_name)
{
  int i, j;
  FILE *pointf;
  char dummy[15];
  int aux[1];
  char aux1;
  int read = 0;


  pointf = fopen (file_name, "r");
  if (pointf == NULL)
    {
      fprintf (stderr, "\nCannot open input file: %s\n", file_name);
      exit (1);
    }

  if (fscanf (pointf, "%s\n", dummy) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }

  //*******************************
  do
    {
      aux1 = getc (pointf);
      //fflush(stdout);
      if (aux1 == '#')
	{
	  do
	    {
	      aux1 = getc (pointf);
	    }
	  while ((aux1 != '\n') && (aux1 != EOF));
	}
      else
	{
	  ungetc (aux1, pointf);
	  read = 1;
	}
    }
  while (read == 0);



  if (fscanf (pointf, " %d %d", xsize, ysize) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }
  //********************************

  if (fscanf (pointf, "%d\n", aux) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }

  for (i = 0; i < *ysize; i++)
    {
      for (j = 0; j < *xsize; j++)
	{
	  pelimg[i][j] = (unsigned char) fgetc (pointf);
	}
    }

  fclose (pointf);		/* closes file */
}

/**
 *
 * <p> Copy the image to memory </p>
 *
 * @param pelimg vector where the image will be saved
 * @param ysize image vertical dimension
 * @param xsize image horizontal dimensio
 * @param file_name file name of the image that will be coded
 */
void
v_read_file_pgm (int *pelimg, int *ysize, int *xsize, char *file_name)
{
  int i, j;
  FILE *pointf;
  char dummy[15];
  int aux[1];
  char aux1;
  int read = 0;


  pointf = fopen (file_name, "r");
  if (pointf == NULL)
    {
      fprintf (stderr, "\nCannot open input file: %s\n", file_name);
      exit (1);
    }

  if (fscanf (pointf, "%s\n", dummy) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }

  //*******************************
  do
    {
      aux1 = getc (pointf);
      //fflush(stdout);
      if (aux1 == '#')
	{
	  do
	    {
	      aux1 = getc (pointf);
	    }
	  while ((aux1 != '\n') && (aux1 != EOF));
	}
      else
	{
	  ungetc (aux1, pointf);
	  read = 1;
	}
    }
  while (read == 0);



  if (fscanf (pointf, " %d %d", xsize, ysize) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }
  //********************************

  if (fscanf (pointf, "%d\n", aux) == EOF)
    {
      printf ("\nEOF ERROR\n");

    }

  for (i = 0; i < *ysize; i++)
    {
      for (j = 0; j < *xsize; j++)
	{
	  pelimg[i * (*xsize) + j] = (unsigned char) fgetc (pointf);
	}
    }

  fclose (pointf);		/* closes file */
}

/**
 * <p> Allocates memory for a matrix of variables of type int. </p>
 *
 * @param nr number of rows
 * @param nc number of columns
 * @return a pointer to a int matrix (int **)
 */
int **
int_matrix (int nr, int nc)
{
  int i;
  int **m;

  m = (int **) malloc ((unsigned) (nr) * sizeof (int *));
  if (!m)
    {
      printf ("int_matrix() - allocation failure 1 \n");
      exit (1);
    }

  for (i = 0; i < nr; i++)
    {
      m[i] = (int *) malloc ((unsigned) (nc) * sizeof (int));
      if (!m[i])
	{
	  printf ("int_matrix() - allocation failure 2 \n");
	  exit (1);
	}
    }

  return m;
}

/**
 * <p> Allocates memory for a matrix of variables of type float. </p>
 *
 * @param nr number of rows
 * @param nc number of columns
 * @return a pointer to a int matrix (float **)
 */
float **
floatmatrix (int nr, int nc)
{
  int i;
  float **m;

  m = (float **) malloc (nr * nc * sizeof (float *));
  if (!m)
    {
      printf ("floatmatrix() - allocation failure 1 \n");
      exit (1);
    }

  for (i = 0; i < nr; i++)
    {
      m[i] = (float *) calloc (nc, sizeof (float));
      if (!m[i])
	{
	  printf ("floatmatrix() - allocation failure 2 \n");
	  exit (1);
	}
    }

  return m;
}

/**
 * <p> Allocates memory for a vector of variables of type int. </p>
 *
 * @param nr number of rows
 * @param nc number of columns
 * @return a pointer to a int vector (int *)
 */
int *
int_vector (int nr, int nc)
{
  int *v;

  v = (int *) malloc ((unsigned) (nr * nc) * sizeof (int *));
  if (!v)
    {
      printf ("int_vector() - allocation failure 1 \n");
      exit (1);
    }

  return v;
}

/**
 * <p> This function will sorte all the pixels of each block to be side by side. </p>
 *
 * @param v_pgm vector with the image
 * @param v_sort_pgm vector with the image blocks sorted (all pixels of the block side by side)
 * @param block_x horizontal size of the block
 * @param block_y vertical size of the block
 * @param xsize horizontal size of the image
 * @param ysize vertical size of the image
 */
void
sort_pgm_blocks (int *v_pgm, int *v_sort_pgm, int block_x, int block_y,
		 int xsize, int ysize)
{
  int count = 0, i, i1, j, j1;

  for (i = 0; i < ysize; i += block_y)
    {
      for (j = 0; j < xsize; j += block_x)
	{

	  //Para todos os blocos
	  for (i1 = 0; i1 < block_y; i1++)
	    {
	      for (j1 = 0; j1 < block_x; j1++)
		{

		  v_sort_pgm[count] = v_pgm[(i + i1) * xsize + (j + j1)];
		  count++;
		}
	    }
	}
    }

}

/**
 * <p> Calculate the Peak Signal Noise Ratio </p>
 *
 * @param origblk
 * @param cmpblk
 * @param nline
 * @param npixel
 * @return
 */
double
calculate_psnr (int **origblk, int **cmpblk, int nline, int npixel)
{
  int i, j;
  double psnr;

  psnr = 0.;
  for (j = 0; j < nline; j++)
    for (i = 0; i < npixel; i++)
      psnr +=
	((double) *(origblk[j] + i) -
	 (double) *(cmpblk[j] + i)) * ((double) *(origblk[j] + i) -
				       (double) *(cmpblk[j] + i));

  if (psnr == 0.0)
    psnr = 0.000000000000000000000000000000000000000000000000000001;
  psnr =
    10.0 * log10 ((double) (RANGEY * RANGEY) /
		  (psnr / (double) (nline * npixel)));

  return psnr;
}

/**
 * <p> Mean Squared Error </p>
 *
 * @param origblk
 * @param cmpblk
 * @param nline
 * @param npixel
 * @return
 */
double
calculate_mse (int **origblk, int **cmpblk, int nline, int npixel)
{
  int i, j;
  long cnt = 0;
  double mse;

  mse = 0.;
  for (j = 0; j < nline; j++)
    for (i = 0; i < npixel; i++)
      {
	mse +=
	  ((double) *(origblk[j] + i) -
	   (double) *(cmpblk[j] + i)) * ((double) *(origblk[j] + i) -
					 (double) *(cmpblk[j] + i));
	cnt++;
      }

  return (mse / cnt);
}

/**
 * <p> Output a bit </p>
 *
 * @param bit
 * @param output_file
 * @param buffer
 * @param bits_to_go
 * @param bits_count
 */
void
output_bit (int bit, FILE * output_file, int *buffer, int *bits_to_go,
	    long *bits_count)
{
  *buffer >>= 1;		/* Put bit in top of buffer  */
  if (bit)
    *buffer |= 0x80;
  *bits_to_go -= 1;
  *bits_count += 1;
  if (*bits_to_go == 0)		/* Output buffer if it is     */
    {				/* now full                                        */
      putc (*buffer, output_file);
      /*       printf("  %x\n", buffer); */
      *bits_to_go = 8;
    }
}

/**
 * <p> Flush out the last bits </p>
 *
 * @param output_file
 * @param buffer
 * @param bits_to_go
 */
void
done_outputing_bits (FILE * output_file, int *buffer, int *bits_to_go)
{
  putc (*buffer >> *bits_to_go, output_file);
  //fprintf(stderr," Total Bits: %d bits\n",Bits_Count);
}

/**
 * <p> Write output to luminance file - PGM format </p>
 *
 * @param im_matrix
 * @param nline
 * @param npixel
 * @param filename
 */
void
write_f_pgm (int **im_matrix, int nline, int npixel, char *filename)
{
  int i;
  int pointfo;
  char header_pgm[20];
  int npixel_orig, nline_orig;
  //ssize_t bytes_written;

  npixel_orig = npixel;
  nline_orig = nline;

  if ((pointfo = creat (filename, PERMS)) == -1)	/* creates file */
    {
      fprintf (stderr, "Write_f_PGM - cannot create PGM file -> %s",
	       filename);
      exit (1);
    }

  /* Defines the header of file - type pgm */
  *(header_pgm + 0) = 'P';
  *(header_pgm + 1) = '5';
  *(header_pgm + 2) = '\n';

  i = 3;
  if (npixel >= 1000)
    {
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

  if (nline >= 1000)
    {
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

  if (write (pointfo, (char *) header_pgm, i) == -1)
    {
      printf ("something went wrong");
    }





  unsigned char **image_tmp;
  int j;
  image_tmp = ucmatrix (0, nline, 0, npixel);

  for (i = 0; i < nline; i++)
    for (j = 0; j < npixel; j++)
      image_tmp[i][j] = (unsigned char) im_matrix[i][j];

  for (i = 0; i < (nline_orig); i++)
    {
      if (write (pointfo, image_tmp[i], npixel_orig) == -1)
	{
	  printf ("something went wrong on %d", i);
	}
    }
  close (pointfo);		/* closes file */
}

/**
 * <p> Allocates memory for a matrix of variables of type unsigned char </p>
 *
 * @param nrl number of last row
 * @param nrh number of first row
 * @param ncl number of last column
 * @param nch number of first column
 * @return a pointer to a unsigned char matrix (unsigned char **)
 */
unsigned char **
ucmatrix (int nrl, int nrh, int ncl, int nch)
{
  int i;
  unsigned char **m;

  m =
    (unsigned char **) malloc ((unsigned) (nrh - nrl + 1) * sizeof (char *));
  if (!m)
    {
      printf ("ucmatrix() - allocation failure 1 \n");
      exit (1);
    }
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
    {
      m[i] =
	(unsigned char *) malloc ((unsigned) (nch - ncl + 1) * sizeof (char));
      if (!m[i])
	{
	  printf ("ucmatrix() - allocation failure 2 \n");
	  exit (1);
	}
      m[i] -= ncl;
    }
  return m;
}

////////////////////////////////////////////////////////////////////////////////
///   Copyright (C) 2008 by Nelson Carreira Francisco                        ///
///   eng.nelsito@gmail.com                                                  ///
///                                                                          ///
///   This program is free software; you can redistribute it and/or modify   ///
///   it under the terms of the GNU General Public License as published by   ///
///   the Free Software Foundation; either version 2 of the License, or      ///
///   (at your option) any later version.                                    ///
///                                                                          ///
///   This program is distributed in the hope that it will be useful,        ///
///   but WITHOUT ANY WARRANTY; without even the implied warranty of         ///
///   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ///
///   GNU General Public License for more details.                           ///
///                                                                          ///
///   You should have received a copy of the GNU General Public License      ///
///   along with this program; if not, write to the                          ///
///   Free Software Foundation, Inc.,                                        ///
///   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.              ///
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
///   Implementacaoo de um codificador de imagens baseado em                 ///
///   Quantificacao vectorial                                                ///
///   Nelson Carreira Francisco                                              ///
////////////////////////////////////////////////////////////////////////////////


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
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "params.h"


////////////////////////////////////////////////////////////////////////////////
///                         CONSTANTS ERROR                                  ///
////////////////////////////////////////////////////////////////////////////////
#define ERROR_INVALID_PARAMETERS 1
#define ERROR_ALLOCATE_MEMORY 2
#define ERROR_OPEN_FILE
#define ERROR_INVALID_PARAMETERS 1
