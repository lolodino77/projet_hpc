#include <stdlib.h>
#include <stdio.h>
#include <time.h>		/* chronometrage */
#include <string.h>		/* pour memset */
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#include "rasterfile.h"

char info[] = "\
Usage:\n\
      mandel dimx dimy xmin ymin xmax ymax prof\n\
\n\
      dimx,dimy : dimensions de l'image a generer\n\
      xmin,ymin,xmax,ymax : domaine a calculer dans le plan complexe\n\
      prof : nombre maximale d'iteration\n\
\n\
Quelques exemples d'execution\n\
      mandel 800 800 0.35 0.355 0.353 0.358 200\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 500\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 300\n\
      mandel 800 800 -1.48478 0.00006 -1.48440 0.00044 100\n\
      mandel 800 800 -1.5 -0.1 -1.3 0.1 10000\n\
";

void print_array(unsigned char* T, int len){
	int i;
	for(i = 0;i < len;i++){
		printf("%d ", T[i]);
	}
	printf("\n\n");
}

double my_gettimeofday()
{
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

/**
 * Convertion entier (4 octets) LINUX en un entier SUN
 * @param i entier à convertir
 * @return entier converti
 */

int swap(int i)
{
	int init = i;
	int conv;
	unsigned char *o, *d;

	o = ((unsigned char *)&init) + 3;
	d = (unsigned char *)&conv;

	*d++ = *o--;
	*d++ = *o--;
	*d++ = *o--;
	*d++ = *o--;

	return conv;
}

/*** 
 * Par Francois-Xavier MOREL (M2 SAR, oct2009): 
 */

unsigned char power_composante(int i, int p)
{
	double iD = i;
	iD /= 255.0;
	iD = pow(iD, p);
	iD *= 255;
	unsigned char o = (unsigned char)iD;
	return o;
}

unsigned char cos_composante(int i, double freq)
{
	double iD = i;
	iD = cos(iD / 255.0 * 2 * M_PI * freq);
	iD += 1;
	iD *= 128;
	unsigned char o = (unsigned char)iD;
	return o;
}

/*** 
 * Choix du coloriage : definir une (et une seule) des constantes
 * ci-dessous :  
 */
//#define ORIGINAL_COLOR
#define COS_COLOR

#ifdef ORIGINAL_COLOR
#define COMPOSANTE_ROUGE(i)    ((i)/2)
#define COMPOSANTE_VERT(i)     ((i)%190)
#define COMPOSANTE_BLEU(i)     (((i)%120) * 2)
#endif				/* #ifdef ORIGINAL_COLOR */
#ifdef COS_COLOR
#define COMPOSANTE_ROUGE(i)    cos_composante(i,13.0)
#define COMPOSANTE_VERT(i)     cos_composante(i,5.0)
#define COMPOSANTE_BLEU(i)     cos_composante(i+10,7.0)
#endif				/* #ifdef COS_COLOR */

unsigned char xy2color(double a, double b, int prof)
{
	double x = 0;
	double y = 0;
	int i;
	for (i = 0; i < prof; i++) {
		/* garder la valeur précédente de x qui va etre ecrase */
		double temp = x;
		/* nouvelles valeurs de x et y */
		double x2 = x * x;
		double y2 = y * y;
		x = x2 - y2 + a;
		y = 2 * temp * y + b;
		if (x2 + y2 > 4.0)
			break;
	}
	return (i == prof) ? 255 : (int)((i % 255));
}

void sauver_rasterfile(char *nom, int largeur, int hauteur, unsigned char *p)
{
	FILE *fd;
	struct rasterfile file;
	int i;
	unsigned char o;

	if ((fd = fopen(nom, "w")) == NULL) {
		printf("erreur dans la creation du fichier %s \n", nom);
		exit(1);
	}

	file.ras_magic = swap(RAS_MAGIC);
	file.ras_width = swap(largeur);	/* largeur en pixels de l'image */
	file.ras_height = swap(hauteur);	/* hauteur en pixels de l'image */
	file.ras_depth = swap(8);	/* profondeur de chaque pixel (1, 8 ou 24 )   */
	file.ras_length = swap(largeur * hauteur);	/* taille de l'image en nb de bytes          */
	file.ras_type = swap(RT_STANDARD);	/* type de fichier */
	file.ras_maptype = swap(RMT_EQUAL_RGB);
	file.ras_maplength = swap(256 * 3);

	fwrite(&file, sizeof(struct rasterfile), 1, fd);

	/* Palette de couleurs : composante rouge */
	i = 256;
	while (i--) {
		o = COMPOSANTE_ROUGE(i);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	/* Palette de couleurs : composante verte */
	i = 256;
	while (i--) {
		o = COMPOSANTE_VERT(i);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	/* Palette de couleurs : composante bleu */
	i = 256;
	while (i--) {
		o = COMPOSANTE_BLEU(i);
		fwrite(&o, sizeof(unsigned char), 1, fd);
	}

	// pour verifier l'ordre des lignes dans l'image : 
	//fwrite( p, largeur*hauteur/3, sizeof(unsigned char), fd);

	// pour voir la couleur du '0' :
	// memset (p, 0, largeur*hauteur);

	fwrite(p, largeur * hauteur, sizeof(unsigned char), fd);
	fclose(fd);
}


//h_div + 1 au lieu de h_div
int main(int argc, char* argv[]){
	int my_rank;
	int p;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	enum tagType {INDICE, TRAITEMENT, STOP};
	unsigned char* Image;
	unsigned char* ImageK;
	MPI_Status status;
	int tag = 0;
	int dest;
	int source;

	int w, h;
	int prof;
	int i, j;
	int k;
	int a, b;
	int numBloc = 0;
	int bTmp = 0;
	int idTmp;
	double x, y;
	double xmin, ymin;
	double xmax, ymax;
	double xinc, yinc;
	unsigned char *ima;
	double debut, fin;
	int tagFin;

	debut = my_gettimeofday();

	xmin = -2;
	ymin = -2;
	xmax = 2;
	ymax = 2;
	w = h = 800;
	prof = 10000;
	int h_div = 80;
	xinc = (xmax - xmin) / (w - 1);
	yinc = (ymax - ymin) / (h - 1);

	Image = (unsigned char*)malloc(sizeof(unsigned char)*w*h);
	ImageK = (unsigned char*)malloc(sizeof(unsigned char)*w*h_div);

	if(my_rank == 0){
		/* Remplissage de la partie de l'image assignee au processus 0 */
		y = ymin;

		/* Premier tour des ordres envoyes aux esclaves */
		for(i = 1;i < p;i++){
			dest = i;
			MPI_Send(&numBloc, 1, MPI_INT, dest, INDICE, MPI_COMM_WORLD);
			numBloc += 1;			
		}

		/* Envoi des ordres */
		while(numBloc != h/h_div){

			////////BUG LA
			MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
			MPI_Recv(ImageK, sizeof(unsigned char)*w*h_div, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
			idTmp = status.MPI_SOURCE;
			dest = status.MPI_SOURCE;
			MPI_Send(&numBloc, 1, MPI_INT, dest, INDICE, MPI_COMM_WORLD);
			numBloc += 1;
			
			/* Remplissage grace au travail d'un sous-tableau */
			for(i = 0;i<w*h_div;i++){
				if(i == (w*h_div-1)){
				}
				Image[bTmp*w*h_div + i] = ImageK[i];
			}
		}

		/* Reception des derniers travaux des esclaves */
		for(i = 1;i < p;i++){
			MPI_Recv(&bTmp, 1, MPI_INT, MPI_ANY_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);	
			MPI_Recv(ImageK, sizeof(unsigned char)*w*h_div, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, TRAITEMENT, MPI_COMM_WORLD, &status);
			idTmp = status.MPI_SOURCE;
			for(j= 0;j<w*(h_div);j++){
				(Image + bTmp*w*h_div)[j] = ImageK[j];
			}	
			MPI_Send(&xmin, 1, MPI_INT, idTmp, STOP, MPI_COMM_WORLD);			
		}
		sauver_rasterfile("mandel.ras", w, h, Image);
	}
	else{
		while(1){
			MPI_Recv(&numBloc, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			tagFin = status.MPI_TAG;
			if(tagFin == STOP){
				break;
			}

			/* Calcul */
			y = ymin + numBloc*(h_div)*yinc;
			for (i = 0; i < (h_div); i++){
				x = xmin;
				for (j = 0; j < w; j++){
					ImageK[j + i * w] = xy2color(x, y, prof);
					//max ImageK[w+w*(h_div)]
					x += xinc;
				}
				y += yinc;
			}
			dest = 0;
			bTmp = numBloc;
			MPI_Send(&bTmp, 1, MPI_INT, dest, TRAITEMENT, MPI_COMM_WORLD);
			MPI_Send(ImageK, sizeof(unsigned char)*w*((h_div)), MPI_UNSIGNED_CHAR, dest, TRAITEMENT, MPI_COMM_WORLD);	
		}
	}

	/* Affichage de la sortie */
	fin = my_gettimeofday();
	fprintf(stderr, "Temps total de calcul du processeur %d : %g sec\n", my_rank, fin - debut);
	//printf("my_rank = %d, fini\n",my_rank);

	MPI_Finalize();
}
