#include <pthread.h>
#include <ctime>

void synchronize(int total_threads);
double get_full_time();

typedef struct
{
  //input
  int n;
  double *matrix;
  double *reverse;
  double *max_el_list;
  int *max_idx_list;
  int *idxs;
  int thread_idx;
  int total_threads;
  //output
  int *algo_error;
  double thread_time;
} ARGS;
