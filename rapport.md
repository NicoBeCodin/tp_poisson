# Exercice 3 : Référence et utilisation de BLAS/LAPACK

## 1. En C, comment doit-on déclarer et allouer une matrice pour utiliser BLAS et LAPACK ?

Dans BLAS et LAPACK, les matrices sont stockées en mémoire sous la forme de tableaux 1D (linéaires) 
en utilisant le format Colonne-Major (Column-Major Order).

### Déclaration et allocation :
```c
int m = 3, n = 3; // Dimensions de la matrice (3x3 par exemple)
double *A = (double *)malloc(m * n * sizeof(double)); // Allocation dynamique d'une matrice 3x3
```
Une matrice 3x3 sera stocké dans la mémoire: [a11​,a21​,a31​,a12​,a22​,a32​,a13​,a23​,a33​]

```c
A[0] = 1.0; // a11
A[1] = 2.0; // a21
A[2] = 3.0; // a31
A[3] = 4.0; // a12
//et ainsi de suite...
free(A)//Libération de la mémoire
```

## Quelle est la signification de la constante LAPACK_COL_MAJOR ?
La constante LAPACK_COL_MAJOR est utilisée pour indiquer que les matrices sont stockées en ordre colonne-major (Column-Major Order),
le format natif de Fortran, le langage pour lequel BLAS/LAPACK a été conçu. Les éléments d'une matrice sont donc stockés colonne par
colonne, ce qui diffère du C, où le stockage est fait ligne par ligne.


## 3. À quoi correspond la dimension principale (leading dimension) généralement notée ld ?

La leading dimension correspond au nombre d'éléments entre deux élements consécutifs de la même
colonn dans la mémoire. Avec la constante LAPACK_COL_MAJOR, la leading dimension correspond au nombre de lignes
de la matrice. La leading dimesnion est donc nécéssaire dans l'appel d'une fonction LAPACK ou BLAS pour assurer un 
accès correct aux données.

## 4. Que fait la fonction dgbmv ? Quelle méthode implémente t-elle ?
La fonction dgbmv effectue la multiplication d'une matrice stockée sous forme de bandes
par un vecteur. Elle implémente le produit matrice-vecteur.

## 5. Que fait la fonction dgbtrf ? Quelle méthode implémente t-elle ?
La fonction dgbtrf réalise la factorisation LU d'une matrice stockée sous forme de bandes.
Elle implémente la factorisation LU et peut donc résoudre des systèmes linéaires ou calculer le determinant d'une matrice

## 6. Que fait la fonction dgbtrs ? Quelle méthode implémente t-elle ?
La fonction dgbtrs résout un système linéaire avec une matrice stockée sous forme de bandes en utilisant
le résulatat de la factorisation LU efefctué par dgbtrf. Elle implémente la méthode de substitution arrière
après la factorisation LU.

## 7. Que fait la fonction dgbsv ? Quelle méthode implémente t-elle ?
La fonction dgbsv résout un système linéaire avec une matrice stockée sous forme de bandes en faisant 
la factorisation LU et la substitution arrière en une étape. Elle implémente la méthode de factorisation LU avec 
substitution arrière.

## 8. Comment calculer la norme du résidu relatif avec des appels BLAS ?

On peut utiliser la fonction dnrm2 pour calculer la norme euclidienne d'un vecteur. On peut l'appliquer au vecteur du résidu après avoir
effectué une opération matrice-vecteur avec dgbmv et obtenu la solution du système linéaire.


