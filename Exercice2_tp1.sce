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
