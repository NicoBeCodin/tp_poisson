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
le résulatat de la factorisation LU effectué par dgbtrf. Elle implémente la méthode de substitution arrière
après la factorisation LU.

## 7. Que fait la fonction dgbsv ? Quelle méthode implémente t-elle ?
La fonction dgbsv résout un système linéaire avec une matrice stockée sous forme de bandes en faisant 
la factorisation LU et la substitution arrière en une étape. Elle implémente la méthode de factorisation LU avec 
substitution arrière.

## 8. Comment calculer la norme du résidu relatif avec des appels BLAS ?

On peut utiliser la fonction dnrm2 pour calculer la norme euclidienne d'un vecteur. On peut l'appliquer au vecteur du résidu après avoir
effectué une opération matrice-vecteur avec dgbmv et obtenu la solution du système linéaire.







### Complexité algorithmique, en espace, et qualité numérique

| **Fonction**                        | **Complexité temporelle** | **Complexité en espace** | **Qualité numérique** | **Commentaire**                                                                                           |
|-------------------------------------|---------------------------|---------------------------|------------------------|-----------------------------------------------------------------------------------------------------------|
| `set_GB_operator_colMajor_poisson1D` | \( O(la) \)               | \( O(la \cdot kv) \)      | Bonne                 | Allocation pour stocker une matrice en bande (\( kv \) dépend de la largeur de bande, souvent constante). |
| `set_GB_operator_colMajor_poisson1D_Id` | \( O(la) \)            | \( O(la \cdot kv) \)      | Bonne                 | Similaire à la fonction précédente, qualité numérique robuste car il n'y a pas de calcul sensible.        |
| `set_dense_RHS_DBC_1D`              | \( O(la) \)               | \( O(la) \)               | Bonne                 | Gère les conditions aux bords avec des initialisations simples, donc stable.                             |
| `set_analytical_solution_DBC_1D`    | \( O(la) \)               | \( O(la) \)               | Bonne                 | Stable, mais dépend de la précision des valeurs initiales pour les conditions aux bords.                  |
| `set_grid_points_1D`                | \( O(la) \)               | \( O(la) \)               | Bonne                 | Crée une grille uniforme, très stable numériquement.                                                      |
| `relative_forward_error`            | \( O(la) \)               | \( O(1) \)                | Dépendante            | Calcul des erreurs sensibles aux petites perturbations dans les données, dépend de la précision des données. |
| `indexABCol`                        | \( O(1) \)                | \( O(1) \)                | Bonne                 | Simple calcul d'indice, insensible aux erreurs numériques.                                                 |
| `dgbtrftridiag`                     | \( O(1) \)                | \( O(1) \)                | N/A                   | Fonction placeholder sans calcul réel.                                                                    |
| `eig_poisson1D`                     | \( O(1) \)                | \( O(1) \)                | Bonne                 | Calcule un scalaire en combinant des fonctions stables.                                                   |
| `eigmax_poisson1D`                  | \( O(1) \)                | \( O(1) \)                | Bonne                 | Qualité robuste car les calculs sont basés sur des opérations mathématiques simples.                       |
| `eigmin_poisson1D`                  | \( O(1) \)                | \( O(1) \)                | Bonne                 | Idem que pour `eigmax_poisson1D`.                                                                         |
| `richardson_alpha_opt`              | \( O(1) \)                | \( O(1) \)                | Bonne                 | Aucun problème de stabilité vu la simplicité des calculs.                                                 |
| `richardson_alpha`                  | \( O(k \cdot la) \)       | \( O(la) \)               | Moyenne à Bonne       | Qualité dépendante du choix de \( \alpha \) et de la convergence des itérations.                          |
| `richardson_MB`                     | \( O(k \cdot la) \)       | \( O(la) \)               | Moyenne à Bonne       | Sensible aux erreurs d'approximation dues à la matrice de bande modifiée.                                  |
| `extract_MB_jacobi_tridiag`         | \( O(la) \)               | \( O(la) \)               | Bonne                 | Robuste car il s'agit d'une extraction directe sans calcul complexe.                                       |
| `extract_MB_gauss_seidel`           | \( O(la) \)               | \( O(la) \)               | Bonne                 | Idem que pour `extract_MB_jacobi_tridiag`, très stable.                                                   |
| `set_CSR_poisson1D`                 | \( O(N) \)                | \( O(3N - 2) \)           | Bonne                 | Génération efficace d'une matrice creuse au format CSR, sans risque significatif d'erreurs numériques.     |
| `set_CSC_poisson1D`                 | \( O(N) \)                | \( O(3N - 2) \)           | Bonne                 | Génération efficace d'une matrice creuse au format CSC, également stable numériquement.                   |
| `dcsrmv`                            | \( O(N) \)                | \( O(N) \)                | Bonne                 | Produit matrice-vecteur au format CSR, stable numériquement grâce à des opérations simples.                |
| `dcscmv`                            | \( O(N) \)                | \( O(N) \)                | Bonne                 | Produit matrice-vecteur au format CSC, stable numériquement et efficace pour les matrices creuses.         |


---

### Justifications des choix

1. **Complexité temporelle**:
   - Les fonctions qui parcourent une grille ou une matrice ont \( O(la) \).
   - Les fonctions de calcul élémentaire (comme des opérations scalaires) ont \( O(1) \).
   - Les méthodes directes de la bibliotheques lapacke ont \( O(n^3) \)
   - Méthodes itératives : \( O(k \cdot n) \) où \( k \) est le nombre d'itérations.
   - Les fonctions itératives dépendent du nombre d'itérations \( k \), d'où \( O(k \cdot la) \).
   - `dcsrmv` et `dcscmv` réalisent des produits matrice-vecteur, chaque entrée non nulle de la matrice étant multipliée et accumulée, d'où \( O(N) \).

2. **Complexité en espace**:
   - Les fonctions allouant de la mémoire pour des matrices (par exemple, stockées en bande ou dense) ont une complexité \( O(la) \) ou \( O(la \cdot kv) \).
   - Les fonctions qui n'allouent pas de grandes structures (comme `indexABCol`) ont \( O(1) \).
   - Les formats CSR et CSC stockent uniquement les valeurs non nulles (\( 3N - 2 \) pour une matrice de Poisson 1D), avec des indices supplémentaires.

3. **Qualité numérique**:
   - Les fonctions de configuration de matrices ou vecteurs (e.g., `set_` fonctions) sont stables, car elles manipulent des données déterministes.
   - Les calculs itératifs (comme `richardson_alpha`) sont sensibles aux erreurs d'arrondi et aux mauvais choix de paramètres (\( \alpha \)).
   - Les erreurs relatives (`relative_forward_error`) dépendent de la précision des données et de l'échelle des erreurs, donc leur qualité peut varier.


---
---

## Résultats et Convergence

### Graphes de convergence

Les images suivantes montrent la convergence des méthodes en fonction des itérations :

#### Méthode de Richardson (alpha)
![Convergence Richardson (alpha)](resvec_iter_0.png)

- ~30 cycles cpu avec les paramètres par défaut

#### Méthode de Jacobi
![Convergence Jacobi](resvec_iter_1.png)

- ~100 cyles cpu avec les paramètres par défaut

#### Méthode de Gauss-Seidel
![Convergence Gauss-Seidel](resvec_iter_2.png)

- ~80 cycles cpu avec les paramètres par défaut
---

### Comparaison des méthodes utilisées

Trois méthodes itératives ont été utilisées :
- **Richardson (alpha)** : Méthode simple basée sur un paramètre optimal.
- **Jacobi** : Méthode basée sur une matrice diagonale préconditionnée.
- **Gauss-Seidel** : Méthode basée sur une précondition \( (D-E) \), converge plus rapidement que Jacobi pour les matrices diagonales dominantes.


### Méthodes directes et itératives

| **Méthode**        | **Avantages**                                                                 | **Inconvénients**                                                                                         |
|---------------------|------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------|
| Méthodes Directes  | - Précision élevée<br>- Convient pour les matrices denses<br>- Résolution exacte | - Complexité temporelle élevée \(O(n^3)\)<br>- Consomme beaucoup de mémoire                                                                   |
| Méthodes Itératives | - Scalabilité pour les matrices creuses<br>- Complexité temporelle \(O(k \cdot n)\) | - Convergence dépendante de la condition de la matrice<br>- Sensibles aux erreurs d'arrondi pour certaines configurations |



## Conclusion

1. Les méthodes itératives comme Jacobi et Gauss-Seidel sont bien adaptées aux grandes matrices creuses.
2. Gauss-Seidel est généralement plus performant que Jacobi en termes de convergence.
3. La méthode de Richardson peut être utilisée comme point de départ pour comprendre les méthodes itératives.



