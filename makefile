# simple make file
SOURCES=CO2-Ar_PES_GP_Symm.f90 PES_GP_Symm.f90
PRODUCT=CO2-Ar_PES.out


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	gfortran -o $(PRODUCT) $(SOURCES)
