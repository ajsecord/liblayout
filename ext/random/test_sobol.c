#include <stdio.h>
#include <stdlib.h>
#include "random.h"

int main(int argc, char* argv[]) {
	long i, j;
	long num_random;
	int num_group;
	float* x;

	if (argc < 3) { 
		fprintf(stderr, "Usage: %s num_random num_group\n", argv[0]);
		exit(1);
	} else {
		num_random = atol(argv[1]);
		num_group = atol(argv[2]);
	}

	x = malloc(num_group * sizeof(float));
	if (!x) {
		fprintf(stderr, "Could not allocate array of %li bytes\n", 
				num_group * sizeof(float));
		exit(1);
	}

	/* Initialise the Sobol generator */
	num_group *= -1;
	rng_sobol(&num_group, x);
	num_group *= -1;

	/* Generate the points, num_group at a time */
	for (i = 0; i < num_random / num_group; i++) {
		rng_sobol(&num_group, x);
		for (j = 0; j < num_group; j++) {
			fprintf(stdout, "%G ", x[j]);
		}
		fprintf(stdout, "\n");
	}

	free(x);

	return 0;
}
