#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

void assign_vector_to_matrix(gsl_matrix* mat, gsl_vector* v, int index, int isRow) {
	if (isRow == 1) {
		for (uint32_t i=0; i<mat->size2; i++) {
			gsl_matrix_set(mat, index, i, v->data[i]);
		}
	} else {
		for (uint32_t i=0; i<mat->size1; i++) {
			gsl_matrix_set(mat, i, index, v->data[i]);
		}
	}
}

void assign_value_to_vec(gsl_vector* v, double* val) {
	for (uint32_t i=0; i<v->size; i++) {
		v->data[i] = val[i];
	}
}

double norm2(gsl_vector* v) {
	double sum = 0;
	for (uint32_t i=0; i<v->size; i++) {
		sum += v->data[i] * v->data[i];
	}

	return sqrt(sum);
}
double norm2Diff(gsl_vector* v1, gsl_vector* v2) {
	gsl_vector* newV = gsl_vector_alloc(v1->size);
	for (uint32_t i=0; i<newV->size; i++) {
		newV->data[i] = v1->data[i] - v2->data[i];
	}
	double result = norm2(newV);
	gsl_vector_free(newV);

	return result;
}
gsl_matrix* invert_a_matrix(gsl_matrix* matrix) {
	size_t size = matrix->size1;

	gsl_permutation *p = gsl_permutation_alloc(size);

	int s;

	//Compute the LU decoposition of this matrix
	gsl_linalg_LU_decomp(matrix, p, &s);
	
	//Compute the inverse of the LU decomposition
	gsl_matrix* inv = gsl_matrix_alloc(size, size);
	gsl_linalg_LU_invert(matrix, p, inv);

	gsl_permutation_free(p);

	return inv;
}

gsl_matrix* multiplication_matrix(gsl_matrix* leftMatrix, gsl_matrix* rightMatrix) {
	gsl_matrix* result = gsl_matrix_alloc(leftMatrix->size1, rightMatrix->size2);

	for (uint32_t i=0; i<leftMatrix->size1; i++) {
		for (uint32_t j=0; j<rightMatrix->size2; j++) {
			double sum = 0;
			for (uint32_t k=0; k<rightMatrix->size1; k++) {
				sum += gsl_matrix_get(leftMatrix, i, k) * gsl_matrix_get(rightMatrix, k, j);
			}	
			gsl_matrix_set(result, i, j, sum);
		}
	}
	return result;
}

gsl_vector* multiplication_matrix_and_vec(gsl_matrix* mat, gsl_vector* v) {
	gsl_vector* result = gsl_vector_alloc(mat->size1);

	for (uint32_t i=0; i<mat->size1; i++) {
			double sum = 0;
			for (uint32_t k=0; k<v->size; k++) {
				sum += gsl_matrix_get(mat, i, k) * v->data[k];
			}	
			result->data[i] = sum;
	}
	return result;
}

gsl_matrix* p_invert_a_matrix(gsl_matrix* matrix) {
	size_t size = matrix->size2;

	gsl_matrix* t_matrix = gsl_matrix_alloc(matrix->size2, matrix->size1);
	gsl_matrix* multiplication = gsl_matrix_alloc(size, size);
	gsl_matrix* result_matrix = gsl_matrix_alloc(matrix->size2, matrix->size1);

	// make transpose matrix
	for (uint32_t i=0; i<matrix->size1; i++) {
		for (uint32_t j=0; j<matrix->size2; j++) {
			gsl_matrix_set(t_matrix, j, i, gsl_matrix_get(matrix, i, j));
		}
	}

	// show the transpose matrix
	//for (uint32_t i=0; i<t_matrix->size1; i++) {
	//	for (uint32_t j=0; j<t_matrix->size2; j++) {
	//		printf("%f ", gsl_matrix_get(t_matrix, i, j));
	//	}
	//	printf("\n");
	//}

	// execute the matrix multiplication
	for (uint32_t i=0; i<t_matrix->size1; i++) {
		double sum = 0;
		for (uint32_t j=0; j<matrix->size2; j++) {
			for (uint32_t k=0; k<matrix->size1; k++) {
				sum += gsl_matrix_get(t_matrix, i, k) * gsl_matrix_get(matrix, k, j);
			}
			gsl_matrix_set(multiplication, i, j, sum);
			sum = 0;
		}
	}

	// show the multipilication result
	//for (uint32_t i=0; i<size; i++) {
	//	for (uint32_t j=0; j<size; j++) {
	//		printf("%f ", gsl_matrix_get(multiplication, i, j));
	//	}
	//	printf("\n");
	//}

	// get the invert of the multiplication matrix
	gsl_matrix* inv_multiplication = invert_a_matrix(multiplication);

	// execute multiplication
	for (uint32_t i=0; i<inv_multiplication->size1; i++) {
		double sum = 0;
		for (uint32_t j=0; j<t_matrix->size2; j++) {
			for (uint32_t k=0; k<t_matrix->size1; k++) {
				sum += gsl_matrix_get(inv_multiplication, i, k) * gsl_matrix_get(t_matrix, k, j);
			}
			gsl_matrix_set(result_matrix, i, j, sum);
			sum = 0;
		}
	}
	
	// show the p_invert matrix result
	//for (uint32_t i=0; i<result_matrix->size1; i++) {
	//	for (uint32_t j=0; j<result_matrix->size2; j++) {
	//		printf("%f ", gsl_matrix_get(result_matrix, i, j));
	//	}
	//	printf("\n");
	//}


	return result_matrix;
}


void show_matrix(gsl_matrix* mat) {
	for (uint32_t i=0; i<mat->size1; i++) {
		for (uint32_t j=0; j<mat->size2; j++) {
			printf("%10f ", gsl_matrix_get(mat, i, j));
		}
		printf("\n");
	}
}

void show_vector(gsl_vector* v) {
	for (uint32_t i=0; i<v->size; i++) {
		printf("%10f ", v->data[i]);
	}
	printf("\n");
}

void initial_vector(gsl_vector* v, double val) {
	for (uint32_t i=0; i<v->size; i++) {
		v->data[i] = val;
	}
}
