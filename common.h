#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
void assign_vector_to_matrix(gsl_matrix* mat, gsl_vector* v, int index, int isRow) ;
void assign_value_to_vec(gsl_vector* v, double* val) ;
double norm2(gsl_vector* v) ;
double norm2Diff(gsl_vector* v1, gsl_vector* v2) ;
gsl_matrix* invert_a_matrix(gsl_matrix* matrix) ;
gsl_matrix* multiplication_matrix(gsl_matrix* leftMatrix, gsl_matrix* rightMatrix) ;
gsl_vector* multiplication_matrix_and_vec(gsl_matrix* mat, gsl_vector* v) ;
gsl_matrix* p_invert_a_matrix(gsl_matrix* matrix) ;
void initial_vector(gsl_vector* v, double val) ;

void show_matrix(gsl_matrix* mat) ;
void show_vector(gsl_vector* v) ;

#endif
