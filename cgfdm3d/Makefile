#================================================================
#   Copyright (C) 2021 Sangfor Ltd. All rights reserved.
#   
#   File Name：Makefile
#   Author: Wenqiang Wang
#   Created Time:2021-10-30
#   Discription:
#
#================================================================

GPU_CUDA := ON 


#FREE_SURFACE := ON
#PML := ON 


CCHOME := /usr
CUDAHOME := /public/software/cuda-10.0
MPIHOME := /public/software/openmpi-4.1.1-cuda



CC := $(CCHOME)/bin/gcc

#General Compiler
ifdef GPU_CUDA
GC := $(CUDAHOME)/bin/nvcc -Xptxas=-v -arch=sm_70 -rdc=true -maxrregcount=128 
else
GC := $(CCHOME)/bin/g++
endif


LIBS := -L$(CUDAHOME)/lib64 -lcudart -lcublas
INCS := -I$(CUDAHOME)/include 

LIBS += -L$(MPIHOME)/lib -lmpi
INCS += -I$(MPIHOME)/include 


SRCDIR := ./src
OBJDIR := ./obj
BINDIR := ./bin


CFLAGS := -c -O2 
LFLAGS := -O2

GCFLAGS := 

ifdef GPU_CUDA
#LFLAGS += -Xptxas=-v 

#LFLAGS += -arch=sm_70 -rdc=true -Xptxas=-v 
#GCFLAGS += --fmad=false 
GCFLAGS += -x cu 
endif

vpath

vpath % $(SRCDIR)
vpath % $(OBJDIR)
vpath % $(BINDIR)


DFLAGS_LIST := GPU_CUDA \
			   FREE_SURFACE PML


DFLAGS := $(foreach flag,$(DFLAGS_LIST),$(if $($(flag)),-D$(flag)))


OBJS := cjson.o printInfo.o readParams.o \
		init_gpu.o cpu_Malloc.o init_grid.o init_MPI.o \
		run.o \
		coord.o medium.o calc_CFL.o\
		contravariant.o \
		wave_deriv.o wave_rk.o freeSurface.o \
		init_pml_para.o pml_deriv.o pml_freeSurface.o \
		propagate.o \
		pml_rk.o \
		source.o\
		data_io.o \
		MPI_send_recv.o MPI_send_recv_jac.o \
		PGV.o \
		main.o 

OBJS := $(addprefix $(OBJDIR)/,$(OBJS))


$(BINDIR)/main : $(OBJS)
	$(GC) $(LFLAGS) $(LIBS) $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(GC) $(CFLAGS) $(DFLAGS) $(GCFLAGS) $(INCS)  $^ -o $@


$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	-rm $(OBJDIR)/* -rf
	-rm $(BINDIR)/* -rf
	-rm output/* -rf


