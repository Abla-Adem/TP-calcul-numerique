CC=gcc
Fprog=-std=c99 -O3 -o
LIBS=-llapacke -lcblas -lm
c="100"
l="5"
n="5.123"
type="N"
alpha="1"
beta="1"
all:exec exec_lu exec_g
exec_lu:Lu.o
	$(CC)  -o exe Lu.o $(OPTC)  $(LIBS)
main.o: Lu.c
	$(CC) $(Fprog) Lu.o -c Lu.c -Wall -O
exec_g:m.o
	$(CC)  -o ex L.o $(OPTC)  $(LIBS)
m.o: Gauss_jacobi.c
	$(CC) $(Fprog) L.o -c Gauss_jacobi.c -Wall -O

exec:mai.o
	$(CC)  -o exec Exercice4_tp3.o $(OPTC)  $(LIBS)
mai.o: Exercice4_tp3.c
	$(CC) $(Fprog) Exercice4_tp3.o -c Exercice4_tp3.c -Wall -O
run:
	./exec $(c) $(l) $(n) $(type) $(alpha) $(beta)
run_lu:
	./exe $(c) $(l) $(n) $(type) $(alpha) $(beta)
rung:
	./ex $(c) $(alpha)
clean:
	rm *.o

