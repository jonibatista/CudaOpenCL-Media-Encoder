#define G_BlocksPerGrid   65536
#define G_ThreadsPerBlock 1024

__kernel void
encoding_pgm (const int num_codewords, const int pgm_block_size,
	      const int num_blocks, const int num_threads,
	      __global int *dev_dict, __global int *dev_pgm,
	      __global int *dev_pgm_coded)
{
  __local float cache_err[G_ThreadsPerBlock];
  __local int cache_idx[G_ThreadsPerBlock];


  int i, idx_dict, idx_block;
  int tid = get_global_id (0);
  int jump = G_BlocksPerGrid * G_ThreadsPerBlock;
  float temp = 0.0;

  if (get_local_id (0) == 0)
    {
      for (i = 0; i < G_ThreadsPerBlock; i++)
	{
	  cache_err[i] = FLT_MAX;
	  cache_idx[i] = 0;
	}
    }

  barrier (CLK_LOCAL_MEM_FENCE);

// make the dictionary stride  
  while (tid < (G_BlocksPerGrid * num_codewords))

    {

      i = 0;
      temp = 123.0;
      idx_dict = get_local_id (0) * pgm_block_size;
      idx_block = 0;

      while (i < pgm_block_size)
	{
	  idx_block = (get_group_id (0) * pgm_block_size) + i;
	  temp +=
	    ((dev_dict[idx_dict + i] -
	      dev_pgm[idx_block]) * (dev_dict[idx_dict + i] -
				     dev_pgm[idx_block]));
	  i++;
	}

      if (cache_err[get_local_id (0)] > temp)
	{
	  cache_err[get_local_id (0)] = temp;
	  cache_idx[get_local_id (0)] = idx_dict/pgm_block_size;
	}
      tid += jump;
    }

  barrier (CLK_LOCAL_MEM_FENCE);

  if (get_local_id (0) == 0)
    {
      float aux = FLT_MAX;
 
      for (i = 0; i < G_ThreadsPerBlock; i++)
	{
	  if (cache_err[i] < aux)
	    {
	      aux = cache_err[i];
	      dev_pgm_coded[get_group_id (0)] = cache_idx[i];
	    }
	}
    }

  barrier (CLK_GLOBAL_MEM_FENCE);
}
