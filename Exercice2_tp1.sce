function [erreur_avant,erreur_arriere,cnd,x]=erreur(n)
a=rand(n,n)
xex=rand(n,1)
b=a*xex
x=a\b
erreur_avant=norm(x-xex)/norm(xex)
erreur_arriere=norm(b-a*x)/norm(b)
cnd=cond(a)
x=cnd*erreur_arriere
endfunction
//Definir n
[erreur_avant,erreur_arriere,cnd,err]=erreur(n)
/*
n=10000
 erreur_avant  = 

   4.301778323D-11
 erreur_arriere  = 

   6.766716024D-15
 cnd  = 

   2.378639409D+06
 err  = 

   1.609557740D-08
   */
