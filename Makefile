TARGET_EXEC := Hyph-m

FC:=gfortran

BIN_DIR := ./bin
BUILD_DIR := ./build
MOD_DIR := ./build/mod
SRC_DIR := ./src
SRC_EXEC := $(TARGET_EXEC).F90

# For Ubuntu+brew
FCFLAGS:=-O3 -I/home/linuxbrew/.linuxbrew/Cellar/fortran-stdlib/0.8.1/include/fortran_stdlib/GNU-15.2.0 -lfortran_stdlib
LDFLAGS:=-L/home/linuxbrew/.linuxbrew/lib -L/home/linuxbrew/.linuxbrew/Cellar/fortran-stdlib/0.8.1/lib -lfortran_stdlib -lfortran_stdlib_stats -lfortran_stdlib_system -lfortran_stdlib_stringlist -lfortran_stdlib_specialmatrices -lfortran_stdlib_quadrature -lfortran_stdlib_logger -lfortran_stdlib_linalg_iterative -lfortran_stdlib_io -lfortran_stdlib_hashmaps -lfortran_stdlib_ansi -lfortran_stdlib_sparse -lfortran_stdlib_specialfunctions -lfortran_stdlib_selection -lfortran_stdlib_math -lfortran_stdlib_linalg -lfortran_stdlib_sorting -lfortran_stdlib_bitsets -lfortran_stdlib_lapack_extended -lfortran_stdlib_lapack -lfortran_stdlib_strings -lfortran_stdlib_intrinsics -lfortran_stdlib_blas -lfortran_stdlib_linalg_core -lfortran_stdlib_hash -lfortran_stdlib_constants -lfortran_stdlib_core -lfortran_stdlib_array
LDFLAGS_LINALG:=-L/home/linuxbrew/.linuxbrew/lib -llapack
FCFLAGS_LINALG:=
# For MacOs+brew
# FCFLAGS:=-O3 -I/opt/homebrew/Cellar/fortran-stdlib/0.8.1/include -I/opt/homebrew/Cellar/fortran-stdlib/0.8.1/include/fortran_stdlib/GNU-15.2.0
# LDFLAGS:=-L/opt/homebrew/Cellar/fortran-stdlib/0.8.1/lib -lfortran_stdlib -lfortran_stdlib_stats -lfortran_stdlib_system -lfortran_stdlib_stringlist -lfortran_stdlib_specialmatrices -lfortran_stdlib_quadrature -lfortran_stdlib_logger -lfortran_stdlib_linalg_iterative -lfortran_stdlib_io -lfortran_stdlib_hashmaps -lfortran_stdlib_ansi -lfortran_stdlib_sparse -lfortran_stdlib_specialfunctions -lfortran_stdlib_selection -lfortran_stdlib_math -lfortran_stdlib_linalg -lfortran_stdlib_sorting -lfortran_stdlib_bitsets -lfortran_stdlib_lapack_extended -lfortran_stdlib_lapack -lfortran_stdlib_strings -ldl -lm -lfortran_stdlib_intrinsics -lfortran_stdlib_blas -lfortran_stdlib_linalg_core -lfortran_stdlib_hash -lfortran_stdlib_constants -lfortran_stdlib_core -lfortran_stdlib_array
# FCFLAGS_LINALG:=-I/opt/homebrew/Cellar/openblas/0.3.34/include -I/opt/homebrew/opt/libomp/include -Xpreprocessor -fopenmp
# LDFLAGS_LINALG:=-L/opt/homebrew/Cellar/openblas/0.3.34/lib -lopenblas

LDFLAGS:=$(LDFLAGS) $(LDFLAGS_LINALG)
EXTRA_FCFLAGS:=-g -ffree-line-length-none -Wall -Wextra -Wpedantic\
               -Wcharacter-truncation\
               -Wfrontend-loop-interchange\
               -Wrealloc-lhs -Wrealloc-lhs-all\
               -funroll-loops -fcheck=bounds\
               -finit-integer=-666 -finit-real=nan\
               -fgcse-lm -fgcse-sm -ftree-vectorize\
               -ffpe-trap=invalid,zero,overflow # -Warray-temporaries -Wconversion-extra
FCFLAGS:=$(FCFLAGS) $(FCFLAGS_LINALG) $(EXTRA_FCFLAGS)
FCFLAGS+= -J $(MOD_DIR)

# Find all the F90 files we want to compile
SRCS := $(shell find $(SRC_DIR) -name '*.F90')

# Object files
OBJS := $(subst $(SRC_DIR),$(BUILD_DIR),$(SRCS)) # subst src dir with build dir
OBJS := $(OBJS:.F90=.o) # change extensions
OBJ_EXEC := $(BUILD_DIR)/$(SRC_EXEC)
OBJ_EXEC := $(OBJ_EXEC:.F90=.o)
$(info $(OBJS))
$(info $(OBJ_EXEC))

# Compile exec
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(FC) $(OBJS) -o $@ $(LDFLAGS)

# Build objects
$(BUILD_DIR)/commons.o: $(SRC_DIR)/commons.F90
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/readInput.o: $(SRC_DIR)/readInput.F90
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/lapack_wrapper.o: $(SRC_DIR)/lapack_wrapper.F90 $(BUILD_DIR)/commons.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/dynmat.o: $(SRC_DIR)/dynmat.F90 $(BUILD_DIR)/commons.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/phonons.o: $(SRC_DIR)/phonons.F90 $(BUILD_DIR)/commons.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/solve_saparov.o: $(SRC_DIR)/solve_saparov.F90 $(BUILD_DIR)/commons.o $(BUILD_DIR)/lapack_wrapper.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/solve_phmag.o: $(SRC_DIR)/solve_phmag.F90 $(BUILD_DIR)/commons.o $(BUILD_DIR)/lapack_wrapper.o $(BUILD_DIR)/dynmat.o $(BUILD_DIR)/phonons.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

$(BUILD_DIR)/main.o: $(SRC_DIR)/main.F90 $(BUILD_DIR)/commons.o $(BUILD_DIR)/lapack_wrapper.o $(BUILD_DIR)/dynmat.o $(BUILD_DIR)/phonons.o $(BUILD_DIR)/solve_phmag.o
	$(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
	$(shell test -d $(MOD_DIR) || mkdir -p $(MOD_DIR))
	$(FC) -c $< -o $@ $(FCFLAGS)

# Install
install: $(BUILD_DIR)/$(TARGET_EXEC)
	$(shell test -d $(BIN_DIR) || mkdir -p $(BIN_DIR))
	$(shell ln -s ../$(BUILD_DIR)/$(TARGET_EXEC) $(BIN_DIR)/)

# Clean
.PHONY: clean
clean:
	if [ -e $(BUILD_DIR) ]; then\
	  rm -r $(BUILD_DIR) &> /dev/null;\
	fi
	if [ -e $(BIN_DIR) ]; then\
	  rm -r $(BIN_DIR) &> /dev/null;\
	fi
