make:
	touch s.txt
	#rm s*

	mpicc -o a.out project3.c -fopenmp -lm
	squeue
	sbatch -N8 jobscript.sh
