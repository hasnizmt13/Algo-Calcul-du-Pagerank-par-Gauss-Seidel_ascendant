#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Définition des constantes
#define fichier_resultat "resultats/resultats_gauss.txt"

// Variables globales
int quantite_memoire_allouee = 0;
int quantite_memoire_liberee = 0;
int nbr_iteration_convergence = 0;

// Structures de données

// Représente un élément dans une liste chaînée, stockant la provenance et la probabilité de lien
struct element {
    int prov;
    double proba;
    struct element *suiv;
};
typedef struct element ELEMENT;

// Représente une ligne dans une matrice, contenant une liste d'éléments
struct ligne {
    int num;
    ELEMENT *elem;
};
typedef struct ligne LIGNE;

// Représente une matrice avec un tableau de lignes
struct matrice {
    int nbr_elem_non_nul;
    int nbr_lignes;
    LIGNE *ligne;
};

// Initialisation des listes chaînées de chaque ligne de la matrice
struct matrice init_matrice(struct matrice matrice) {
    for(int i = 0; i < matrice.nbr_lignes; i++) {
        matrice.ligne[i].elem = NULL;
    }
    return matrice;
}

// Ajoute un nouvel élément en tête de la liste chaînée
ELEMENT *nouvel_element(ELEMENT *elem, int prov, double proba) {
    ELEMENT *nouvel_elem = malloc(sizeof(ELEMENT));
    quantite_memoire_allouee += sizeof(ELEMENT);
    nouvel_elem->prov = prov;
    nouvel_elem->proba = proba;
    nouvel_elem->suiv = elem;
    return nouvel_elem;
}

// Libération de la mémoire allouée pour une liste chaînée
ELEMENT *liberation_liste(ELEMENT *elem) {
    if(elem != NULL) {
        elem->suiv = liberation_liste(elem->suiv);
        free(elem);
        quantite_memoire_liberee += sizeof(ELEMENT);
    }
    return NULL;
}

// Libération de la mémoire allouée pour une matrice
void liberation_matrice(struct matrice matrice) {
    for(int i = 0; i < matrice.nbr_lignes; i++) {
        matrice.ligne[i].elem = liberation_liste(matrice.ligne[i].elem);
    }
    quantite_memoire_liberee += sizeof(LIGNE) * matrice.nbr_lignes;
    free(matrice.ligne);
}

// Lecture des données d'un fichier et création de la matrice
struct matrice lecture(char const *nom_fichier, int stanford) {
    struct matrice matrice;
    int degre = 0;
    int dest = 0;
    double proba = 0.0;

    FILE *fichier = fopen(nom_fichier, "r");
    if(fichier == NULL) {
        printf("Erreur : Le fichier renseigné est introuvable\n");
        exit(0);
    }
    printf("\nLecture du fichier...\n");

    fscanf(fichier, "%d\n", &matrice.nbr_elem_non_nul);
    fscanf(fichier, "%d\n", &matrice.nbr_lignes);

    matrice.ligne = malloc(sizeof(LIGNE) * matrice.nbr_lignes);
    quantite_memoire_allouee += sizeof(LIGNE) * matrice.nbr_lignes;

    matrice = init_matrice(matrice);

    for(int i = 0; i < matrice.nbr_lignes; i++) {
        fscanf(fichier, "%d %d ", &matrice.ligne[i].num, &degre);
        if(stanford) matrice.ligne[i].num--;

        for(int j = 0; j < degre; j++) {
            if(stanford) {
                fscanf(fichier, "%d %lf ", &dest, &proba);
                matrice.ligne[dest-1].elem = nouvel_element(matrice.ligne[dest-1].elem, i, proba);
            } else {
                fscanf(fichier, "%lf %d ", &proba, &dest);
                matrice.ligne[dest].elem = nouvel_element(matrice.ligne[dest].elem, i, proba);
            }
        }
        printf("\r");
        printf("%d %%", i * 100 / matrice.nbr_lignes + 1);
    }
    printf("\n\n");

    fclose(fichier);
    return matrice;
}

// Écriture des résultats dans un fichier
void ecriture_resultat(double *pin, int taille) {
    FILE *fichier = fopen(fichier_resultat, "w+");
    if(fichier == NULL) {
        printf("Erreur : Il est impossible d'écrire dans le fichier %s\n", fichier_resultat);
    } else {
        fprintf(fichier, "(%lf", pin[0]);
        for(int i = 1; i < taille; i++) {
            fprintf(fichier, ", %lf", pin[i]);
        }
        fprintf(fichier, ")\n");
        fclose(fichier);
    }
}

// Calcul du PageRank en utilisant la méthode de Gauss-Seidel Ascendant
double *pagerank_gauss_seidel_ascendant(struct matrice matrice) {
    double eps = 0.000001;
    double abs = 1.0;
    double tmp;
    double alpha = 0.85;
    double *pio, *pin;
    ELEMENT *elem_parcours;
    double somme_renormalisation = 0;

    pio = malloc(sizeof(double) * matrice.nbr_lignes);
    quantite_memoire_allouee += sizeof(double) * matrice.nbr_lignes;

    pin = malloc(sizeof(double) * matrice.nbr_lignes);
    quantite_memoire_allouee += sizeof(double) * matrice.nbr_lignes;

    for (int i = 0; i < matrice.nbr_lignes; ++i)
        pio[i] = 1.0 / matrice.nbr_lignes;

    printf("Calcul du pagerank :\n\n");
    while(eps < abs) {
        nbr_iteration_convergence++;
        abs = 0.0;
        somme_renormalisation = 0;

        for (int i = 0; i < matrice.nbr_lignes; ++i) {
            pin[i] = 0.0;
            elem_parcours = matrice.ligne[i].elem;
            while(elem_parcours != NULL) {
                if(elem_parcours->prov < i) 
                    pin[i] += pin[elem_parcours->prov] * elem_parcours->proba;
                else 
                    pin[i] += pio[elem_parcours->prov] * elem_parcours->proba;
                elem_parcours = elem_parcours->suiv;
            }
            pin[i] = (1 - alpha) / matrice.nbr_lignes + alpha * pin[i];
            somme_renormalisation += pin[i];
        }

        for (int i = 0; i < matrice.nbr_lignes; ++i) {
            pin[i] = pin[i] / somme_renormalisation;
            tmp = pin[i] - pio[i];
            if(tmp < 0) tmp = -tmp;
            abs += tmp;
            pio[i] = pin[i];
        }
        printf("Difference entre 2 itérations : %lf\n", abs);
    }

    quantite_memoire_liberee += sizeof(double) * matrice.nbr_lignes;
    free(pio);

    return pin;
}

// Programme principal
int main(int argc, char const *argv[]) {
    int stanford = 0;
    clock_t debut_lecture_t, fin_lecture_t, debut_calcul_t, fin_calcul_t;
    struct matrice matrice;

    if(argc == 3) {
        if(!strcmp(argv[2], "--stanford"))
            stanford = 1;
        else {
            printf("Erreur : L'argument donné n'est pas reconnu.\n");
            exit(0);
        }
    } else if(argc != 2) {
        printf("Erreur : Le nombre d'arguments donnés est incorrect.\n");
        exit(0);
    }

    debut_lecture_t = clock();
    matrice = lecture(argv[1], stanford);
    fin_lecture_t = clock();

    debut_calcul_t = clock();
    double *pin = pagerank_gauss_seidel_ascendant(matrice);
    fin_calcul_t = clock();

    ecriture_resultat(pin, matrice.nbr_lignes);
    quantite_memoire_liberee += sizeof(double) * matrice.nbr_lignes;
    free(pin);

    liberation_matrice(matrice);

    printf("\n*******************************************\n");
    printf("\nNombre d'itérations pour la convergence : %d\n", nbr_iteration_convergence);
    printf("\nTemps de lecture du fichier en CPU ticks : %lu\n", fin_lecture_t - debut_lecture_t);
    printf("Temps de lecture du fichier en ms : %lf\n", (double)(fin_lecture_t - debut_lecture_t) * 1000 / CLOCKS_PER_SEC);
    printf("Temps de calcul du pagerank en CPU ticks : %lu\n", fin_calcul_t - debut_calcul_t);
    printf("Temps de calcul du pagerank en ms : %lf\n", (double)(fin_calcul_t - debut_calcul_t) * 1000 / CLOCKS_PER_SEC);
    printf("\nTemps total de calcul du pagerank en CPU ticks : %lu\n", fin_calcul_t - debut_lecture_t);
    printf("Temps total de calcul du pagerank en ms : %lf\n", (double)(fin_calcul_t - debut_lecture_t) * 1000 / CLOCKS_PER_SEC);
    printf("\nQuantité de mémoire allouée dynamiquement en octets : %d\n", quantite_memoire_allouee);
    printf("Quantité de mémoire libérée dynamiquement en octets : %d\n", quantite_memoire_liberee);

    return 0;
}
