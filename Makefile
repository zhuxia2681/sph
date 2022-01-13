FC=mpiifort -g
CC=mpiicc -g
FLAG= -DPARA -DVERBOSE -DLOG
OBJ=cell.o io.o number.o integral.o numerical.o \
mp.o sli1.o random.o util.o

default: cell io number integral numerical mp sli1\
	random util
	$(FC) $(FLAG) $(INC) scp.F90 $(OBJ)


io: mp
	$(FC) $(FLAG) -c io.F90

mp: 
	$(FC) $(FLAG) -c mp.F90

number: 
	$(FC) $(FLAG) -c number.F90

cell: number
	$(FC) $(FLAG) -c cell.F90

integral: number numerical random
	$(FC) $(FLAG) -c integral.F90

numerical: number
	$(FC) $(FLAG) -c numerical.F90

sli1: number numerical
	$(FC) $(FLAG) -c sli1.F90 

random: random
	$(FC) $(FLAG) -c random.F90 

util: 
	$(FC) $(FLAG) -c util.F90 


clean:
	rm -f *.o *.mod a.out
