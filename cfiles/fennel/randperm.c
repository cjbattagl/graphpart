// randperm.c
// casey battaglino 11/7/2013

#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

// Generate a random permutation of 1:size using the "Knuth shuffle"
int* genRandPerm(int size) {
  int *orderList = (int *) malloc (sizeof (int) * size);
  assert(*orderList);
  srand(time(NULL));
  int randomPick, temp = 0;
  int remainingNumbers = size-1;
  
  // Generate 'identity' permutation
  for (int i = 0; i < size; i++) { orderList[i] = i; }

  // Shuffle
  while (remainingNumbers > 0) {
    randomPick = ((rand() % (remainingNumbers+1)));
    temp = orderList[remainingNumbers];
    orderList[remainingNumbers] = orderList[randomPick];
    orderList[randomPick] = temp;
    remainingNumbers--;
  }
	
  // Test code: assert that all numbers sum to n(n-1)/2
  //long sum = 0;
  //for (int i=0; i < size; i++) { sum += orderList[i]; }
  //assert(sum == ((size-1)*(size))/2);  
  return orderList;
}