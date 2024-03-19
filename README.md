make all : compile tous les fichiers .c
make clean : supprime tous les exécutables

Commande de compilation : gcc -fopenmp -Wall -O3  -o nom_code nom_code.c -lm -mavx2
Commande d’exécution : ./nom_code --size <N> 	(N : puissance de 2) (on utilise 67108864)
Commande changement de nombre de threads : export OMP_NUM_THREADS=n