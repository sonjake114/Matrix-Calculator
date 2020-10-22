#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, double low, double high) {
    srand(42);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    if (rows <= 0 || cols <= 0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    *mat = malloc(sizeof(matrix));
    if (*mat == NULL) {
      free(*mat);
      PyErr_SetString(PyExc_RuntimeError, "Memory error");
      return -1;
    }
    (*mat) -> rows = rows;
    (*mat) -> cols = cols;
    (*mat) -> ref_cnt = 1;
    (*mat) -> parent = NULL;
    int size = rows * cols;

    double * datafile = (double *)calloc(size, sizeof(double));
    if (datafile == NULL){
      free(datafile);
      PyErr_SetString(PyExc_RuntimeError, "Memory error");
      return -1;
    }
    (*mat) -> data  = datafile;
    return 0;

    /* TODO: YOUR CODE HERE */
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`.
 * You should return -1 if either `rows` or `cols` or both are non-positive or if any
 * call to allocate memory in this function fails. Return 0 upon success.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    if (rows <= 0 || cols <= 0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    *mat = malloc(sizeof(matrix));
    if (*mat == NULL) {
      free(*mat);
      PyErr_SetString(PyExc_RuntimeError, "Memory error");
      return -1;
    }
    (*mat) -> rows = rows;
    (*mat) -> cols = cols;
    double * fromdata = from -> data;
    fromdata = fromdata + offset;
    (*mat) -> data = fromdata;
    (*mat) -> parent = from;
    from -> ref_cnt = from -> ref_cnt + 1;
    return 0;
}

/*
 * This function frees the matrix struct pointed to by `mat`. However, you need to make sure that
 * you only free the data if `mat` is not a slice and has no existing slices, or if `mat` is the
 * last existing slice of its parent matrix and its parent matrix has no other references.
 * You cannot assume that mat is not NULL.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat == NULL){
      return;
    }
    if (mat->parent == NULL && mat->ref_cnt == 1){
      free(mat -> data);
      free(mat);
      return;
    }
    if (mat -> parent != NULL && mat->parent->ref_cnt == 2 && mat->ref_cnt == 2){
      free(mat -> data);
      free(mat);
      return;
    }
    if (mat -> parent != NULL && mat->parent->ref_cnt == 1){
      free(mat -> parent -> data);
      free(mat -> parent);
      free(mat);
      return;
    }
    mat->parent->ref_cnt = mat->parent->ref_cnt - 1;
    free(mat);
    return;

}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    int width = mat -> cols;
    double val = mat -> data[width*row + col];
    return val;
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    int width = mat -> cols;
    mat -> data[width*row + col] = val;


}

/*
 * Sets all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
  int mat_col = mat -> cols;
  int mat_row = mat -> rows;
  int mat_size =  mat_col * mat_row;
  for (int i = 0; i<mat_size; i++){
    mat -> data[i] = val;
  }

    /* TODO: YOUR CODE HERE */
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if(result == NULL || mat1 == NULL || mat2 == NULL){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    int result_col = result -> cols;
    int result_row = result -> rows;
    int result_size = result_col * result_row;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }


    int mat1_col = mat1 -> cols;
    int mat1_row = mat1 -> rows;
    if (mat1_col<=0 || mat1_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    int mat2_col = mat2 -> cols;
    int mat2_row = mat2 -> rows;
    if (mat2_col<=0 || mat2_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat1_col != mat2_col || mat2_col!=result_col || mat1_col !=result_col){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat1_row != mat2_row || mat2_row!=result_row || mat1_row !=result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    #pragma omp parallel for
    for(int i = 0; i < result_size / 4 * 4 ; i+=4){
      result -> data[i] = (mat1 -> data[i]) + (mat2 -> data[i]);
      result -> data[i + 1] = (mat1 -> data[i + 1]) + (mat2 -> data[i + 1]);
      result -> data[i + 2] = (mat1 -> data[i + 2]) + (mat2 -> data[i + 2]);
      result -> data[i + 3] = (mat1 -> data[i + 3]) + (mat2 -> data[i + 3]);

    }
    #pragma omp parallel for
    for (int i = result_size / 4 * 4; i < result_size; i++) {
	    result -> data[i] = (mat1 -> data[i]) + (mat2 -> data[i]);
    }

    return 0;

}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if(result == NULL || mat1 == NULL || mat2 == NULL){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_col = result -> cols;
    int result_row = result -> rows;
    int result_size = result_col * result_row;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    int mat1_col = mat1 -> cols;
    int mat1_row = mat1 -> rows;
    if (mat1_col<=0 || mat1_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    int mat2_col = mat2 -> cols;
    int mat2_row = mat2 -> rows;
    if (mat2_col<=0 || mat2_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    if (mat1_col != mat2_col || mat2_col!=result_col || mat1_col !=result_col){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat1_row != mat2_row || mat2_row!=result_row || mat1_row !=result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    #pragma omp parallel for
    for(int i = 0; i < result_size / 4 * 4 ; i+=4){
      result -> data[i] = (mat1 -> data[i]) - (mat2 -> data[i]);
      result -> data[i + 1] = (mat1 -> data[i + 1]) - (mat2 -> data[i + 1]);
      result -> data[i + 2] = (mat1 -> data[i + 2]) - (mat2 -> data[i + 2]);
      result -> data[i + 3] = (mat1 -> data[i + 3]) - (mat2 -> data[i + 3]);

    }
    #pragma omp parallel for
    for (int i = result_size / 4 * 4; i < result_size; i++) {
	    result -> data[i] = (mat1 -> data[i]) - (mat2 -> data[i]);
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */

int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {

  if(result == NULL || mat1 == NULL || mat2 == NULL){
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }
  int mat1_col = mat1 -> cols;
  int mat1_row = mat1 -> rows;
  if (mat1_col<=0 || mat1_row<=0){
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }

  int mat2_col = mat2 -> cols;
  int mat2_row = mat2 -> rows;
  int mat2_size = mat2_col * mat2_row;
  if (mat2_col<=0 || mat2_row<=0){
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }
  if (mat1_col != mat2_row) {
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }
  int result_col = result -> cols;
  int result_row = result -> rows;
  if (result_col<=0 || result_row<=0){
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }

  if (mat1_row != result_row){
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }
  if (mat2_col != result_col) {
    PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
    return -1;
  }


  matrix *fake2 = NULL;
  int as= allocate_matrix(&fake2, mat2_col, mat2_row);
  if(as) {
    deallocate_matrix(fake2);
    PyErr_SetString(PyExc_RuntimeError, "Memory allocation error");
    return -1;
  }
  int x, y, z, u;
  #pragma omp for
  for (x= 0; x < mat2_col; x += 16){
    for (y = 0; y < fake2->cols; y += 16) {
      for (z = x; (z < x + 16) && (z < mat2_col); z++){
        for (u = y; (u < y + 16) && (u < fake2->cols); u++){
          fake2->data[z*(fake2->cols) + u] = mat2->data[u*(mat2->cols) + z];
        }
      }
    }
  }

  #pragma omp parallel
    {
      double sum = 0.0;
      #pragma omp for
      for(int x = 0; x < mat1_row; x++){
        for(int y = 0; y < fake2->rows ; y++) {
          sum = 0;
          for(int z = 0; z < mat1_col; z++) {
            sum = sum + mat1->data[x*mat1->cols+z] * fake2->data[y*fake2->cols+z];

          }
          result->data[x*result->cols+y] = sum;
        }
      }
    }
    deallocate_matrix(fake2);



  //if (mat1->cols == mat1->rows && mat2->cols == mat2->rows && mat1_size == mat2_size)


    //dgemm_avx(n, &fake1, &fake2, result);



  /*#pragma omp parallel
    {
      double sum = 0.0;
      #pragma omp for
      for(int x = 0; x < mat1_row; x++){
        for(int y = 0; y < fake->rows ; y++) {
          sum = 0;
          for(int z = 0; z < mat1_col; z++) {
            sum = sum + get(mat1, x, z) * get(fake, y, z);

          }
          set(result, x, y, sum);
        }
      }
    }*/





  /*for (int i = 0; i < mat2_col; i++){
    for (int j = 0; j < result_col; j++) {
        set(result, i, j, get(mat2, j, i));
    }
  }*/



  return 0;

/*  const int UNROLL = 4;

  for (int i=0;  i<mat2->cols;  i+= UNROLL*4) {
      for (int j=0;  j<mat1->rows;  j++) {
          __m256d c[4];
          for (int x=0;  x<UNROLL;  x++)
              c[x] = _mm256_loadu_pd(result->data+i+x*4+j*result_col );
          for (int k=0;  k<mat1_col;  k++) {
              __m256d b = _mm256_broadcast_sd(mat1->data+k+j*mat1_col);
              for (int x=0;  x<UNROLL;  x++)
                  c[x] = _mm256_add_pd(c[x],
                         _mm256_mul_pd(_mm256_loadu_pd(mat2->data+mat2_col*k+x*4+i), b));
          }
          for (int x=0;  x<UNROLL;  x++)
              _mm256_storeu_pd(result->data+i+x*4+j*result_col, c[x]);
      }
  }
  */
  /*const int UNROLL = 4;
  for (int i=0;  i<mat2_col / (4*UNROLL)*(4*UNROLL);  i+= 4*UNROLL) {
      for (int j=0;  j<mat1_row;  j++) {
          __m256d c[4];
          for (int x=0;  x<UNROLL;  x++)
              c[x] = _mm256_loadu_pd(result->data+j+x*4+i*result_col);
          for (int k=0;  k<mat1_col;  k++) {
              __m256d b = _mm256_broadcast_sd(mat2->data+j+k*mat2_col);
              for (int x=0;  x<UNROLL;  x++)
                  c[x] = _mm256_add_pd(c[x],
                         _mm256_mul_pd(_mm256_loadu_pd(mat1->data+mat1_col*i+x*4+k), b));
          }
          for (int x=0;  x<UNROLL;  x++)
              _mm256_storeu_pd(result->data+j+x*4+i*result_col, c[x]);
      }
  }*/

//dgemm_unroll(2, mat1->data, mat2->data, result->data);

  /*for(int x = 0; x < mat1_row; x+=4){
    for(int y = 0; y < mat2_col; y++) {
      __m256d mat2_fake = {0, 0, 0, 0};
      for(int z = 0; z < mat1_col; z++) {
        mat2_fake = _mm256_add_pd(
          mat2_fake,
          _mm256_mul_pd(
            _mm256_loadu_pd((mat1->data)+z+x*mat1_col),
            _mm256_loadu_pd((mat2->data)+y+z*mat2_col)));
      }
      _mm256_storeu_pd((result->data)+y+x*result_col, mat2_fake);
    }
  }*/

/*for (int i=0;  i < mat1_row / 4*4;  i+=4) {
      for (int j=0;  j<mat2_col;  j++) {
          // c0 = c[i][j]
          __m256d c0 = {0,0,0,0};
          for (int k=0;  k<mat1_col;  k+=1) {
              c0 = _mm256_add_pd(
                      c0,   // c0 += a[i][k] * b[k][j]
                      _mm256_mul_pd(
                          _mm256_loadu_pd((mat1->data)+i*mat1_col+k),
                          _mm256_broadcast_sd((mat2->data)+k*mat2_col+j)));
          }
            _mm256_storeu_pd((result->data)+i*result_col+j, c0);
           // c[i,j] = c0
      }
  }*/

/*double sum = 0.0;
  for(int x = mat1_row/4*4; x < mat1_row; x++){
    for(int y = 0; y < mat2_col; y++) {
      sum = 0;
      for(int z = mat1_col/4*4; z < mat1_col; z++) {
        sum = sum + get(mat1, x, z) * get(mat2, z, y);
        set(result, x, y, sum);
      }
    }
  }*/



  /*void dgemm_unroll(int n, double *A, double *B, double *C) {
      for (int i=0;  i<n;  i+= UNROLL*4) {
          for (int j=0;  j<n;  j++) {
              __m256d c[4];
              for (int x=0;  x<UNROLL;  x++)
                  c[x] = _mm256_loadu_pd(C+i+x*4+j*n);
              for (int k=0;  k<n;  k++) {
                  __m256d b = _mm256_broadcast_sd(B+k+j*n);
                  for (int x=0;  x<UNROLL;  x++)
                      c[x] = _mm256_add_pd(c[x],
                             _mm256_mul_pd(_mm256_loadu_pd(A+n*k+x*4+i), b));
              }
              for (int x=0;  x<UNROLL;  x++)
                  _mm256_storeu_pd(C+i+x*4+j*n, c[x]);
          }
      }
  }*/







  /*for(int i = 0; i < mat1_row; i += UNROLL*4){
    for(int j = 0; j < mat2_col; j++) {
      __m256d cole[4];
      for(int x = 0; x < UNROLL; x++) {
        cole[x] = _mm256_loadu_pd((result->data)+i*result_col+x*4+j);
      for(int k = 0; k < mat1_col; k++){
        __m256d bol = _mm256_broadcast_sd((mat2->data)+k*mat2_col+j);
        for (int x = 0; x < UNROLL; x++)
          cole[x] = _mm256_add_pd(cole[x],
                   _mm256_mul_pd(_mm256_load_pd((mat1->data)+i*mat1_col+x*4+k), bol));

      }
      for (int x = 0; x < UNROLL; x++){
        _mm256_storeu_pd((result->data)+i*result_col+x*4+j, cole[x]);
      }
      }
    }
  }*/





  /*  #pragma omp parallel
    {
      double sum = 0.0;
      #pragma omp for
      for(int x = 0; x < mat1_row; x++){
        for(int y = 0; y < mat2_col; y++) {
          sum = 0;
          for(int z = 0; z < mat1_col; z++) {
            sum = sum + get(mat1, x, z) * get(mat2, z, y);
            set(result, x, y, sum);
          }
        }
      }
    }*/





  }



  /*int result_size = result_row * result_col;
  double sum = 0.0;
  for(int x = 0; x < mat1_row; x+=4){
    for(int y = 0; y < mat2_col; y++) {
      __m256d mat2_fake = {0, 0, 0, 0};
      for(int z = 0; z < mat1_col; z++) {
        mat2_fake = _mm256_add_pd(
          mat2_fake,
          _mm256_mul_pd(
            _mm256_loadu_pd((mat1->data)+z+x*mat1_col),
            _mm256_loadu_pd((mat2->data)+y+z*mat2_col)));
      }
      _mm256_storeu_pd((result->data)+y+x*result_col, mat2_fake);
    }
  }*/









    /*for(int x = 0; x < mat1_row; x+=4){
      for(int y = 0; y < mat2_col; y++) {
        __m256d mat2_fake = {0, 0, 0, 0};
        for(int z = 0; z < mat1_col; z++) {
          mat2_fake = _mm256_add_pd(
            mat2_fake,
            _mm256_mul_pd(
              _mm256_loadu_pd((mat1->data)+x+z*mat1_col),
              _mm256_broadcast_sd((mat2->data)+z+y*mat2_col)));
        }
        _mm256_storeu_pd((result->data)+x+y*result_col, mat2_fake);
      }
    }
    */






/*const int UNROLL = 4;

void dgemm_unroll(int n, double *A, double *B, double *C) {
    for (int i=0;  i<n;  i+= 4*UNROLL) {
        for (int j=0;  j<n;  j++) {
            __m256d c[4];
            c[0] = _mm256_loadu_pd(C+i+j*n);
            c[1] = _mm256_loadu_pd(C+i+4+j*n);
            c[2] = _mm256_loadu_pd(C+i+8+j*n);
            c[3] = _mm256_loadu_pd(C+i+12+j*n);
            for (int k=0;  k<n;  k++) {
                __m256d b = _mm256_broadcast_sd(B+k+j*n);
                c[0] = _mm256_add_pd(c[0], _mm256_mul_pd(_mm256_loadu_pd(A+k*n+i), b));
                c[1] = _mm256_add_pd(c[1], _mm256_mul_pd(_mm256_loadu_pd(A+k*n+4+i), b));
                c[2] = _mm256_add_pd(c[2], _mm256_mul_pd(_mm256_loadu_pd(A+k*n+8+i), b));
                c[3] = _mm256_add_pd(c[3], _mm256_mul_pd(_mm256_loadu_pd(A+k*n+12+i), b));
            }
                _mm256_storeu_pd(C+i+j*n, c[0]);
                _mm256_storeu_pd(C+i+4+j*n, c[1]);
                _mm256_storeu_pd(C+i+8+j*n, c[2]);
                _mm256_storeu_pd(C+i+12+j*n, c[3]);
        }
    }


}
*/

/*double sum = 0.0;
for(int x = 0; x < mat1_row; x++){
  for(int y = 0; y < mat2_col; y++) {
    sum = 0;
    for(int z = 0; z < mat1_col; z++) {
      sum = sum + get(mat1, x, z) * get(mat2, z, y);
      set(result, x, y, sum);
    }
  }
}*/


/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
  if (pow <= 0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    if(result == NULL || mat == NULL ){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_col = result -> cols;
    int result_row = result -> rows;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_size = result_col * result_row;
    int mat_col = mat -> cols;
    int mat_row = mat -> rows;
    if (mat_col<=0 || mat_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_col != mat_row || result_col != result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_col != result_col || mat_row != result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    matrix *fake = NULL;
    allocate_matrix(&fake, result_row, result_col);
    #pragma omp for
    for (int i = 0; i < result_size; i++){
      fake -> data[i] = mat -> data[i];
      result -> data[i] = mat -> data[i];
    }

    int k = 0;
    while (k < pow - 1){
      mul_matrix(result, fake, mat);
      #pragma omp for
      for (int i = 0; i < result_size; i++){
        fake -> data[i] = result -> data[i];
      }
      k++;
    }
     deallocate_matrix(fake);
     return 0;




    /* TODO: YOUR CODE HERE */
    /*if (pow <= 0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    if(result == NULL || mat == NULL ){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_col = result -> cols;
    int result_row = result -> rows;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_size = result_col * result_row;
    int mat_col = mat -> cols;
    int mat_row = mat -> rows;
    if (mat_col<=0 || mat_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int mat_size = mat_col * mat_row;
    if (mat_col != mat_row || result_col != result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_col != result_col || mat_row != result_row){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    matrix *fake = NULL;

    allocate_matrix(&fake, result_row, result_col);


    for (int i = 0; i < result_size; i++){
      fake -> data[i] = mat -> data[i];
    }



    for (int k = 0; k < pow - 1; k++){
      mul_matrix(result, fake, mat);
      for (int i = 0; i < result_size; i++){
        fake -> data[i] = result -> data[i];
      }
    }

    deallocate_matrix(fake);
    return 0;*/

}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if(result == NULL || mat == NULL){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_col = result -> cols;
    int result_row = result -> rows;
    int result_size = result_col * result_row;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }


    int mat_col = mat -> cols;
    int mat_row = mat -> rows;
    if (mat_col<=0 || mat_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }

    if (mat_col != result_col){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_row != result_row ){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    #pragma omp parallel for
    for(int i = 0; i < result_size / 4 * 4 ; i+=4){
      result -> data[i]  = - (double) mat -> data[i];
      result -> data[i + 1]  = - (double) mat -> data[i + 1];
      result -> data[i + 2]  = - (double) mat -> data[i + 2];
      result -> data[i + 3]  = - (double) mat -> data[i + 3];
    }
    #pragma omp parallel for
    for (int i = result_size / 4 * 4; i < result_size; i++) {
	    result -> data[i]  = - (double) mat -> data[i];
    }


    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if(result == NULL || mat == NULL){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_col = result -> cols;
    int result_row = result -> rows;
    if (result_col<=0 || result_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    int result_size = result_col * result_row;
    int mat_col = mat -> cols;
    int mat_row = mat -> rows;
    if (mat_col<=0 || mat_row<=0){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_col != result_col){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    if (mat_row != result_row ){
      PyErr_SetString(PyExc_TypeError, "Incorrect number of elements in list");
      return -1;
    }
    #pragma omp parallel for
    for(int i = 0; i < result_size / 4 * 4 ; i+=4){
      result -> data[i]  = fabs(mat -> data[i]);
      result -> data[i + 1]  = fabs(mat -> data[i + 1]);
      result -> data[i + 2]  = fabs(mat -> data[i + 2]);
      result -> data[i + 3]  = fabs(mat -> data[i + 3]);
    }
    #pragma omp parallel for
    for (int i = result_size / 4 * 4; i < result_size; i++) {
	    result -> data[i]  = fabs(mat -> data[i]);
    }



    return 0;
}
