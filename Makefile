#---------------------------------------Transfer functions for CDM
OPT += -DBBKS_T_K
#OPT += -DEH_T_K
#----------------------------------------whether there is WDM, actually sterile neutrino here, K. Abazajian 2006
#OPT += -DWD_ST                           
#----------------------------------------Choice of lambda nu relation and formation time relation
#OPT += -DLEE                            #Lee et al. 2005
OPT += -DDL                              #Our choice
#----------------------------------------Local type non-gaussianity
OPT += -DLOCAL_NG
#--------------------------------------- Select target computer
SYSTYPE="dlcheng"
#SYSTYPE="ITSC"
#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dalong/Install/gsl/include
GSL_LIBS =  -L/home/dalong/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/gsl-1.14/lib  -Wl,"-R /usr/local/gsl-1.14/lib" 
GSL_INCL =  -I/usr/local/gsl-1.14/include
endif



OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = halo_shape_cal

OBJS   = combine.o formation_time.o init.o local_ng.o main.o norm_power.o norm_prob.o var.o

INCL   = var.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) 


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch
