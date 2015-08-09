#include <stdio.h>
#include <stdlib.h>
#include "random.h"

int main(int argc, char* argv[]) {
	long i, j;
	int num_random;
	int num_group;
	float* x;
	long seed = 0;

	if (argc < 3) { 
		fprintf(stderr, "Usage: %s num_random num_group\n", argv[0]);
		exit(1);
	} else {
		num_random = atol(argv[1]);
		num_group = atol(argv[2]);
	}

	x = malloc(num_random * sizeof(float));
	if (!x) {
		fprintf(stderr, "Could not allocate array of %li bytes\n", 
				num_random * sizeof(float));
		exit(1);
	}

	/* Generate the points */
	rng_jittered_halton(num_group, num_random / num_group, 0.5, &seed, x);

	for (i = 0; i < num_random / num_group; i++) {
		for (j = 0; j < num_group; j++) {
			fprintf(stdout, "%G ", x[i * num_group + j]);
		}
		fprintf(stdout, "\n");
	}

	free(x);

	return 0;
}
