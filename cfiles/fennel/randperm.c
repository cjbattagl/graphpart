// randperm.c
// casey b 11/7/2013

#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include "randperm.h"

// Generate a random permutation of 1:size using the "Knuth shuffle"
int* genRandPerm(int size) {
  int *orderList = (int *) malloc (sizeof (int) * size);
  assert(orderList);
  srand(time(NULL));
  
  // Generate 'identity' permutation
  for (int i = 0; i < size; i++) { orderList[i] = i; }
  
  shuffle_int(orderList, size);
  return orderList;
}

void shuffle_int(int *list, int len) {
	int j;
	int tmp;
	while(len) {
      j = irand(len);
      if (j != len - 1) {
        tmp = list[j];
        list[j] = list[len - 1];
        list[len - 1] = tmp;
      }
    len--;
  }
}

int irand(int n) {
	int r, rand_max = RAND_MAX - (RAND_MAX % n);
	while ((r = rand()) >= rand_max);
	return r / (rand_max / n);
}