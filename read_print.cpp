#include "read_print.h"
#include <algorithm>
#include <ctime>

inline static double f(int k, int n, int i, int j)
{
  switch (k) {
    case 1:
      return (double)n - std::max(i, j) + 1;
    case 2:
      return (double)std::max(i, j);
    case 3:
      return (double)abs(i - j);
    case 4:
      return 1.f / (i + j + 1);
    case 5:
      return (double)i == j;
    default:
      return 0;
  }
}

void fill_matrix(int k, int n, double *matrix)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      matrix[i * n + j] = f(k, n, i, j);
}

void print_matrix(int m, int l, int n, double *matrix)
{
  for (int i = 0; i < std::min(m, l); i ++) {
    for (int j = 0; j < std::min(m, n); j++)
      printf("%10.3e ", matrix[i * n + j]);
    printf("\n");
  }
}

int read_matrix(const char *filename, double *matrix, int n)
{
  FILE *input_file;
  input_file = fopen(filename, "r");
  if(input_file == NULL) {
    return -1;
  }
  for (int i = 0; i < n * n; i++) {
    if (fscanf(input_file, "%lf", matrix + i) != 1) {
      fclose(input_file);
      return -2;
    }
  }
  fclose(input_file);
  return 0;
}
