#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/**********
CONSTANTES
**********/
#define fichier_resultat "resultats/resultats_pagerank.txt"


/*********************
VARIABLES GLOBALES
*********************/

int quantite_memoire_allouee = 0;
int quantite_memoire_liberee = 0;
int nbr_iterations_convergence = 0;


/****************************
DEFINITION DES STRUCTURES
****************************/

// Représentation d'un élément dans une matrice
struct element {
	int dest;
	double proba;
};typedef struct element ELEMENT;

// Représentation d'une ligne dans une matrice
struct ligne {
	int num;
	int degre;
	ELEMENT *elem;
};typedef struct ligne LIGNE;

// Représentation d'une matrice
struct matrice {
	int nbr_elem_non_nul;
	int nbr_lignes;
	LIGNE *ligne;
};


/**************************
LIBERATION DE LA MEMOIRE
**************************/

void liberation_matrice(struct matrice matrice) {
	for(int i=0; i<matrice.nbr_lignes; i++) {
		quantite_memoire_liberee += sizeof(ELEMENT)*matrice.ligne[i].degre;
		free(matrice.ligne[i].elem);
	}
	quantite_memoire_liberee += sizeof(LIGNE)*matrice.nbr_lignes;
	free(matrice.ligne);
}


/*****************************************
FONCTIONS DE LECTURE ET D'ECRITURE
*****************************************/

struct matrice lecture(const char *nom_fichier) {
	struct matrice matrice;

	// Ouverture du fichier
	FILE *fichier = fopen(nom_fichier, "r");
	if(fichier == NULL) {
		printf("Erreur : Impossible de trouver le fichier specifie\n");
		exit(0);
	}

	// Lecture des deux premières lignes du fichier
	fscanf(fichier, "%d\n", &matrice.nbr_elem_non_nul);
	fscanf(fichier, "%d\n", &matrice.nbr_lignes);

	// Allocation mémoire pour les lignes de la matrice
	matrice.ligne = malloc(sizeof(LIGNE)*matrice.nbr_lignes);
	quantite_memoire_allouee += sizeof(LIGNE)*matrice.nbr_lignes;

	// Lecture des données de chaque ligne
	printf("\nLecture du fichier...\n");
	for(int i=0; i<matrice.nbr_lignes; i++) {
		// Lecture des informations de base de chaque ligne
		fscanf(fichier, "%d %d ", &matrice.ligne[i].num, &matrice.ligne[i].degre);
		matrice.ligne[i].num--;

		// Allocation mémoire pour les éléments de chaque ligne
		matrice.ligne[i].elem = malloc(sizeof(ELEMENT)*matrice.ligne[i].degre);
		quantite_memoire_allouee += sizeof(ELEMENT)*matrice.ligne[i].degre;

		// Lecture des éléments de la ligne
		for (int j=0; j<matrice.ligne[i].degre; j++) {
				fscanf(fichier, "%d %lf ", &matrice.ligne[i].elem[j].dest, &matrice.ligne[i].elem[j].proba);
				matrice.ligne[i].elem[j].dest--;
		
		}
		printf("\r%d %%", i*100/matrice.nbr_lignes+1);
	}
	printf("\n\n");

	fclose(fichier);
	return matrice;
}

void ecriture_resultat(LIGNE pin) {
	FILE *fichier = fopen(fichier_resultat, "w+");
	if(fichier == NULL) {
		printf("Erreur : Impossible d ecrire dans le fichier %s\n", fichier_resultat);
	} else {
		fprintf(fichier, "(%lf", pin.elem[0].proba);
		for(int i = 1; i < pin.degre; i++) {
			fprintf(fichier, ", %lf", pin.elem[i].proba);
		}
		fprintf(fichier, ")\n");
		fclose(fichier);
	}
}


/*********************
CALCUL DU PAGERANK
*********************/

LIGNE pagerank(struct matrice matrice) {
	double eps = 0.000001;	// Seuil pour la convergence
	double abs = 1.0; 		// Valeur de la différence absolue entre deux itérations
	double tmp;				// Variable temporaire pour le calcul de la différence
	double alpha = 0.85;	// Facteur de probabilité du surfeur aléatoire

	LIGNE pio, pin;

	// Allocation mémoire pour pio et pin
	pio.num = 0, pio.degre = matrice.nbr_lignes;
	pio.elem = malloc(sizeof(ELEMENT)*pio.degre);
	quantite_memoire_allouee += sizeof(ELEMENT)*pio.degre;

	pin.num = 0, pin.degre = matrice.nbr_lignes;
	pin.elem = malloc(sizeof(ELEMENT)*pin.degre);
	quantite_memoire_allouee += sizeof(ELEMENT)*pin.degre;

	// Initialisation de pio
	for (int i = 0; i < pio.degre; ++i)
		pio.elem[i].proba = 1.0/pio.degre;

	// Boucle de convergence
	printf("Calcul du PageRank :\n\n");
	while(eps < abs) {
		nbr_iterations_convergence++;

		// Initialisation de pin
		for (int i = 0; i < pin.degre; ++i)
			pin.elem[i].proba = 0.0;

		// Calcul des nouvelles valeurs de pin
		for (int i = 0; i < matrice.nbr_lignes; ++i) {
			for(int j = 0; j < matrice.ligne[i].degre; j++) {
				pin.elem[matrice.ligne[i].elem[j].dest].proba += pio.elem[i].proba * matrice.ligne[i].elem[j].proba;
			}
		}

		// Ajout du surfeur aléatoire
		for(int i = 0; i < pin.degre; i++) {
			pin.elem[i].proba = (1-alpha)/pin.degre + alpha*pin.elem[i].proba;
		}

		// Calcul de la différence absolue et mise à jour de pio
		abs = 0.0;
		for(int i = 0; i < pin.degre; i++) {
			tmp = pin.elem[i].proba - pio.elem[i].proba;
			if(tmp < 0) tmp = -tmp;
			abs += tmp;
			pio.elem[i].proba = pin.elem[i].proba;
		}
		printf("Difference entre 2 iterations : %lf\n", abs);
	}

	// Libération mémoire de pio
	quantite_memoire_liberee += sizeof(ELEMENT)*pio.degre;
	free(pio.elem);

	return pin;
}

/******************
PROGRAMME PRINCIPAL
******************/
int main(int argc, char const *argv[]) {
	// Variables pour mesurer le temps de calcul et de lecture
	clock_t debut_lecture_t, fin_lecture_t, debut_calcul_t, fin_calcul_t;
	// Structure pour stocker la matrice
	struct matrice matrice;

	// Vérification des arguments d'entrée
	if(argc != 2) {
		printf("Erreur : Nombre incorrect d arguments.\n");
		exit(0);
	}

	// Lecture du fichier
	debut_lecture_t = clock();
	matrice = lecture(argv[1]);
	fin_lecture_t = clock();

	// Calcul du PageRank
	debut_calcul_t = clock();
	LIGNE pin = pagerank(matrice);
	fin_calcul_t = clock();

	// Écriture des résultats dans un fichier et libération de la mémoire
	ecriture_resultat(pin);
	quantite_memoire_liberee += sizeof(ELEMENT)*matrice.nbr_lignes;
	free(pin.elem);

	// Libération de la mémoire de la matrice
	liberation_matrice(matrice);

	// Affichage des statistiques
	printf("\n*******************************************\n");
	printf("\nNombre d'iterations pour la convergence : %d\n", nbr_iterations_convergence);
	printf("\nTemps de lecture du fichier en CPU ticks : %lu\n", fin_lecture_t - debut_lecture_t);
	printf("Temps de lecture du fichier en ms : %lf\n", (double)(fin_lecture_t - debut_lecture_t)*1000/CLOCKS_PER_SEC);
	printf("Temps de calcul du PageRank en CPU ticks : %lu\n", fin_calcul_t - debut_calcul_t);
	printf("Temps de calcul du PageRank en ms : %lf\n", (double)(fin_calcul_t - debut_calcul_t)*1000/CLOCKS_PER_SEC);
	printf("\nTemps total de calcul du PageRank en CPU ticks : %lu\n", fin_calcul_t - debut_lecture_t);
	printf("Temps total de calcul du PageRank en ms : %lf\n", (double)(fin_calcul_t - debut_lecture_t)*1000/CLOCKS_PER_SEC);
	printf("\nQuantite de memoire allouee dynamiquement en octets : %d\n", quantite_memoire_allouee);
	printf("Quantite de memoire liberee dynamiquement en octets : %d\n", quantite_memoire_liberee);

	return 0;
}
