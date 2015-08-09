#include <stdio.h>
#include <stdlib.h>
#include "random.h"

int main(int argc, char* argv[]) {
	long i;
	long num_random;
	long num_group;
	long seed;

	if (argc < 3) { 
		fprintf(stderr, "Usage: %s num_random num_group\n", argv[0]);
		exit(1);
	} else {
		num_random = atol(argv[1]);
		num_group = atol(argv[2]);
	}

	for (i = 0; i < num_random; i++) {
		if (i % num_group == 0) fprintf(stdout, "\n");
		fprintf(stdout, "%G ", rng_uniform_dev(&seed)); 	
	}

	return 0;
}
