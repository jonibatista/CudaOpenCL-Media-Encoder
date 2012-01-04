__kernel void kermel_quadratic(int num_codewords, int pgm_block_size, int *dev_dict, int *dev_pgm, int *dev_pgm_coded){
  __local float cache[65536];
  __local int cache_idx[G_ThreadsPerBlock];

  int i;
  int tid = get_global_id(0);
  int id = get_local_id(0);
  int jump = G_BlocksPerGrid * G_ThreadsPerBlock;
  int cache_index = id;
  float temp = 0.0;

  if (id == 0)
    {
      for (i = 0; i < G_ThreadsPerBlock; i++)
	cache[i] = 0.0;
	cache_idx[i] = 0;
    }

  barrier();

  // make the dictionary stride  
  while (tid < (G_BlocksPerGrid * num_codewords))

   {

  i = 0;
      temp = 0.0;

      int idx_dict = id * pgm_block_size;
      int idx_block = 0;

      while (i < pgm_block_size)
	{
  
  barrier();
 
idx_block = (get_group_id() * pgm_block_size) + i;
	  temp +=
	    ((dev_dict[idx_dict + i] -
	      dev_pgm[idx_block]) * (dev_dict[idx_dict + i] -
				     dev_pgm[idx_block]));
	  i++;
	}

  barrier();
      
if (cache[cache_index] > temp)
	{
	  cache[cache_index] = temp;
	  cache_idx[cache_index] = idx_dict;
	}
  barrier();
      tid += jump; 
    }

  barrier();

  if (id == 0)
    {
      float aux = FLT_MAX;

      for (i = 0; i < G_ThreadsPerBlock; i++)
	{
//	      printf(" %f =>", cache[i]);
	  if (cache[i] < aux);
	    {
	      // printf(" %d ", i);
	      aux = cache[i];
	      dev_pgm_coded[get_group_id()] = cache_idx[i];
	    }
	}
//               printf(" %d=>%d ", get_group_id(), dev_pgm_coded[get_group_id()]);
    }

  barrier();

}
