## Principe
Appliquer les outils de calcul parallèle vu en cours (OpenMPI et OpenMP) sur à un projet de résolution de système linéaire (avec des matrices creuses) dans le cadre d'un projet scolaire.
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