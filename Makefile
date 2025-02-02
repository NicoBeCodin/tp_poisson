##########################################
# Makefile                               #
# Makefile for the code developed in TP1 #
#                                        #
# T. Dufaud                              #
##########################################
################################
# Variables for this makefile
################################
# 
# -- option to dedicated machine
#
# Rk: You should create a file such as your-machineName.mk
# Follow the example of the machine called "ambre" in the 
# file ambre.mk
#
HOSTNAME?=$(shell hostname)
include $(HOSTNAME).mk

# 
# -- Compiler Option
OPTC=${OPTCLOCAL}

#
# -- Directories
TPDIR=.
TPDIRSRC=$(TPDIR)/src

#
# -- librairies
LIBS=${LIBSLOCAL}

# -- Include directories
INCLATLAS=${INCLUDEBLASLOCAL}
INCL= -I $(TPDIR)/include $(INCLATLAS)


OBJENV= tp_env.o
OBJTP2ITER= lib_poisson1D.o tp_poisson1D_iter.o
OBJTP2DIRECT= lib_poisson1D.o tp_poisson1D_direct.o
#

all: bin/tp_testenv bin/tpPoisson1D_direct bin/tpPoisson1D_iter bin/test_csr_csc

testenv: bin/tp_testenv

tpPoisson1D_iter: bin/tpPoisson1D_iter

tpPoisson1D_direct: bin/tpPoisson1D_direct

tp_env.o: $(TPDIRSRC)/tp_env.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_env.c 

lib_poisson1D.o: $(TPDIRSRC)/lib_poisson1D.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/lib_poisson1D.c 

tp_poisson1D_iter.o: $(TPDIRSRC)/tp_poisson1D_iter.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_poisson1D_iter.c  

tp_poisson1D_direct.o: $(TPDIRSRC)/tp_poisson1D_direct.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/tp_poisson1D_direct.c  

csr_csc.o: $(TPDIRSRC)/csr_csc.c
	$(CC) $(OPTC) -c $(INCL) $(TPDIRSRC)/csr_csc.c


bin/tp_testenv: $(OBJENV) 
	$(CC) -o bin/tp_testenv $(OPTC) $(OBJENV) $(LIBS)

bin/tpPoisson1D_iter: $(OBJTP2ITER)
	$(CC) -o bin/tpPoisson1D_iter $(OPTC) $(OBJTP2ITER) $(LIBS)

bin/tpPoisson1D_direct: $(OBJTP2DIRECT)
	$(CC) -o bin/tpPoisson1D_direct $(OPTC) $(OBJTP2DIRECT) $(LIBS)

bin/test_csr_csc: csr_csc.o
	$(CC) -o bin/test_csr_csc $(OPTC) $(TPDIRSRC)/test_csr_csc.c csr_csc.o $(LIBS)


run_testenv:
	bin/tp_testenv

run_tpPoisson1D_iter:
	bin/tpPoisson1D_iter $(ARGS)

run_tpPoisson1D_direct:
	bin/tpPoisson1D_direct $(ARGS)

run_test_csr_csc:
	bin/test_csr_csc


clean:
	rm *.o bin/*