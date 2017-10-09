# simple make file
SOURCES=PES_GP_Symm.f90 
PRODUCT=CO2-CO_PES.out


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	gfortran -o $(PRODUCT) $(SOURCES)
