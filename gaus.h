#include <cmath>


void gaus_fprop(double *matrix, double *reverse, int n, int *idxs,
                int thread_idx, int total_threads, int *max_idx_list,
                double *max_el_list);
int gaus_bprop(double *matrix, double *reverse, int n,
               int total_threads, int thread_idx);
double error(double *matrix, double *reverse, int n);
void replace(double * pre_reverse, double * or_matrix, int n,
             int *idxs, int total_threads, int thread_idx);
void *gaus_full_algorithm(void *void_args);
