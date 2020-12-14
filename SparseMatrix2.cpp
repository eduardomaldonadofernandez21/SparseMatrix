#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <stdlib.h>
#include <mkl.h>

//1 7 0 0; 0 2 8 0; 5 0 3 9; 0 6 0 4

/*
1. Tomar un ejemplo de matriz utilizada en la práctica 1 y reproducir los
resultados con las funciones MKL.

*/

int showMat(double* mat, int row, int col) {
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			printf("%.4f ", mat[col * i + j]);
		}
		printf("\n");
	}
	printf("\n");
	return 0;
}

int showValue(double* mat, int size) {
	for (int i = 0; i < size; i++) {
		if (mat[i] == 0) {
			printf("\n");
			return i;
		}
		printf("%.0f ", mat[i]);
	}
	printf("\n");
	return 0;
}

int showValue(int* mat, int size) {
	for (int i = 0; i < size; i++) {
		printf("%d ", mat[i]);
	}
	printf("\n");
	return 0;
}

int main(int argc, char* argv[]) {
	const int n = 4;
	const int m = 4;
	double A[16] = { 1,7,0,0, 0,2,8,0, 5,0,3,9, 0,6,0,4 };
	double csr[m*n]{ 0 }; 
	int ja[m*n]{ 0 }; 
	int ia[n+1]{ 0 }; 
	int lda = n;
	int info = 1;
	int job[6] = { 0,0,0,2,16,1};

	printf("Matriz dispersa: \n");
	showMat(A,4,4);

	#pragma warning(suppress : 4996)
	mkl_ddnscsr(job,&n,&m,A,&lda,csr,ja,ia,&info);
	printf("Codificacion CSR mediante mkl_ddnscsr:\n");
	printf("Valores:-------------\n");
	int notZerosLenght = showValue(csr, *(&csr + 1) - csr);
	printf("Indice Columnas:-------------\n");
	showValue(ja, notZerosLenght);
	printf("Filas comprimidas:-------------\n");
	showValue(ia, n+1);
	printf("\n");


	printf("Codificacion COO mediante mkl_dscrcoo:\n");
	double coo[16]{ 0 };
	int jobCOO[6] = { 0,0,0,2,16,3 };
	int nnz = *(&csr + 1) - csr;
	int rowind[16]{ 0 };
	int colind[16]{ 0 };
	#pragma warning(suppress : 4996)
	mkl_dcsrcoo(jobCOO,&n,csr,ja,ia,&nnz,coo, rowind, colind, &info);
	printf("Valores:-------------\n");
	notZerosLenght = showValue(coo, *(&coo + 1) - coo);
	printf("Indice Columnas:-------------\n");
	showValue(colind, notZerosLenght);
	printf("Indice Filas:-------------\n");
	showValue(rowind, notZerosLenght);
	printf("\n");

	//[OPTATIVO 3] Probar conversiones a otros formatos.
	double csc[16]{ 0 };
	int jobCSC[6] = { 0,0,0,2,16,3 };
	int ia1[16]{ 0 };
	int ja1[16]{ 0 };
	#pragma warning(suppress : 4996)
	mkl_dcsrcsc(job,&m,csr,ja,ia,csc,ja1,ia1,&info);
	printf("Codificacion CSC mediante mkl_dcsrcsc:\n");
	printf("Valores:-------------\n");
	notZerosLenght = showValue(csc, *(&csc + 1) - csc);
	printf("Columnas comprimidas:-------------\n");
	showValue(ia1, n + 1);
	printf("Indices Filas:-------------\n");
	showValue(ja1, notZerosLenght);
	printf("\n");
	getchar();
	return 0;
}