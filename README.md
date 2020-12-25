# TP-calcul-numerique
depot des Tp calcul numerique  
dans les fichier Exercice... vous trouverais les implementation des algorithmes  
Dans le rapport les tests on etait fait sur des matrice de taille 3:103  
pour faire les Test veuillez regarder le fichier test ou:  
  -a est le nombre d'interation   
  -n la taille de depart  
  -dans le test de l'exercice 4 on peut faire un choix entre les deux methode 0=lsolve 1=usolve  
  pour la serie 2 on a les fichier:  
    -lib_poisson1D.c avec les modification pour l'exercice 3  
    -Exercice4_TP3 on a plusieur argument dans le makefile pour faire varier le nombre de ligne et de colonne(c,l) la valeur a partir de laquel on genere la valeur aleatoire (a) type pour definir le type d'operation  
    -LU pour l'exercice 5  
    -Exercice6_tp3 pour tester la methode on la fonction TestExo6_tp3_conv(a_n,esp,maxit) qui se trouve dans test.sce ou a_n la taille de la matrice ,esp l'epsilon qu'on veut ateindre ,maxit le nombre max d'iteration pour que la methode de gauss-seidel ou jacobi marche il faut:  
      -si la matrice  est symétrique définie positive1 ;  
      -si A est à diagonale strictement dominante.  
pour jacobi en c dans le makefile on a:  
  -alpha pour connaitre le savoir a partir de quelle valeur on fait varier les valeur de la matrice  
  -c la taille de la matrice  
  -il existe deux fonction test :une qui testes les convergnce et une qui teste les perf  
  -les fonction sont en commentaire sur le main(gauss-jacobi)
    
