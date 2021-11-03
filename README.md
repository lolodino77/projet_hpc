# Projet de Calcul haute performance en première année d'école d'ingénieur

## Problème : résolution de système linéaire avec matrices creuses
Appliquer les outils de calcul parallèle (OpenMPI et OpenMP) vu en cours de calcul haute performance (HPC) à un projet de résolution de système linéaire (avec des matrices creuses) dans le cadre d'un projet scolaire en binôme.
Explication du code dans le PDF des consignes du projet.

## Dossiers
* **mpi** : algorithme avec MPI
* **mpi_et_openmp** : algorithme avec MPI et OpenMP
* **certificats** : certificats qui prouvent l'intégrité de la solution trouvée (qu'on n'a pas "triché" mais bien calculer la solution avec un algorithme)

## Fichiers
* **Makefile** : makefile qui pour compiler le programme
* **checkpoint** : valeurs sauvegardées lors du dernier calcul interrompu (peut être importer par le programme pour éviter de recommencer le calcul)
* **hostile, hostfile2 et hostfile_single** : liste des machines/processus auquel on peut faire appel lors du calcul
* **res** : solution du système obtenue après calcul 
