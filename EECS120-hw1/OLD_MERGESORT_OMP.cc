/**
 *  \file mergesort.cc
 *
 *  \brief Implement your mergesort in this file.
 */
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "sort.hh"
void mergesort(keytype* A, int l, int r, keytype *B, int s);
void merge(keytype *A, int l, int m, int r);
int binarySearch(keytype *A, int l, int r, int n);
void pmerge(keytype *A, int p1, int r1, int p2, int r2, keytype *B, int p3);
int a = 1;

void mySort (int N, keytype* A, keytype* A_out)
{       
	keytype *T = newCopy(N,A);
#pragma omp parallel
#pragma omp single
	mergesort(A, 0, N-1, T, 0);
	printf("\n");
	for (int i = 0; i < N; i++){
		std::cout << A[i] << ' ';
	}	
//	free(T);
}
void mergesort(keytype* A, int l, int r, keytype *B, int s)

{
	//if ((r-1) > 100){
	int n = r - l + 1;
	if (n == 1)
		B[s] = A[l];
	else{
		//keytype *T = newKeys(n);
		int m = (l+r)/2;
		int m2 = m - l + 1;
#pragma omp task shared(T)
		mergesort(A,l,m,T,1);
#pragma omp task shared(T)
		mergesort(A,m+1,r,T,m2+1);
#pragma omp taskwait
		pmerge(T,1,m2,m2+1,n,B,s);			
		//free(T);		
	}
}
/*else{
  int n = r - 1 + 1;
  if (n == 1)
  B[s] = A[l];
  else{
  keytype *T = newKeys(n);
  int m = (l+r)/2;
  int m2 = m - l + 1;
  mergesort(A,l,m,T,1);
  mergesort(A,m+1,r,T,m2+1);
  pmerge(T,1,m2,m2+1,n,B,s);
  }
  }
  }*/
void pmerge(keytype *A, int p1, int r1, int p2, int r2, keytype *B, int p3){
//	printf("p1: %d, r1: %d, p2: %d, r2 %d, p3: %d\n",p1,r1,p2,r2,p3);
	int length1 = r1 - p1 + 1;
	int length2 = r2 - p2 + 1;
	if (length1 < length2){
		std::swap(p1,p2);
		std::swap(r1,r2);
		std::swap(length1,length2);
	} 
	if (length1 == 0){
//		printf("RETURN");
		return;

	}
	int q1 = (p1+r1)/2;
	int q2 = binarySearch(A, p2, r2, A[q1]);
	int q3 = p3 + (q1 - p1) + (q2 - p2);
	B[q3] = A[q1];
#pragma omp task
	pmerge(A, p1, q1-1, p2, q2-1, B, p3);
#pragma omp task 
	pmerge(A, q1+1, r1, q2, r2, B, q3+1);
#pragma omp taskwait
}
void merge(keytype* A, int l, int m, int r)
{
	int a = m - l + 1; /*size of left array */
	int b = r - m; /*size of right array*/
	int L[a];
	int R[b];

	for (int i = 0; i < a; i++)
		L[i] = A[l + i]; /*fills the left array*/




	for (int j = 0; j < b; j++)
		R[j] = A[m + 1+ j]; /*fills the right array*/
	int z = l;
	int j = 0;
	int i = 0;

	if (j == -5)
		printf("how?");

	while ((i < a) && (j < b)) /*merges the arrays*/

	{
		if (L[i] >  R[j])
		{
			A[z] = R[j];
			j++;
		}
		else
		{
			A[z] = L[i];
			i++;
		}
		z++;
	}
	if (i == -1)
		printf("error");
	while (i < a)
	{
		A[z] = L[i];
		i++; /* increases i each time*/
		z++;
	}

	while (j < b)
	{
		A[z] = R[j];
		j++;
		z++;
	}
}
int binarySearch(keytype *A, int l, int r, int n){
	long low = l;
	long high = std::max(l,r+1);
	while (low < high){
		long mid = (low+high)/2;
		if (n <= A[mid]) high = mid;
		else low = mid+1;
	}
	//printf("binary search returns %ld\n",high);
	return high;
} 
