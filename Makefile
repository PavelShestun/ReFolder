FC = gfortran
FFLAGS = -O2
LDFLAGS = -llapack -lblas
SRC = src/types.f90 src/io.f90 src/internal_coordinates.f90 src/potential.f90 \
      src/matrices.f90 src/eigenproblem.f90 src/modify_structure.f90 src/metrics.f90 \
      src/main.f90
OBJ = $(SRC:.f90=.o)
EXEC = refolder
TEST_SRC = test_alignment.f90
TEST_OBJ = $(TEST_SRC:.f90=.o)
TEST_EXEC = test_alignment

all: $(EXEC)

$(EXEC): $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS)

src/types.o: src/types.f90
	$(FC) $(FFLAGS) -c $< -o $@

src/io.o: src/io.f90 src/types.o
	$(FC) $(FFLAGS) -c $< -o $@

src/internal_coordinates.o: src/internal_coordinates.f90 src/types.o
	$(FC) $(FFLAGS) -c $< -o $@

src/potential.o: src/potential.f90 src/types.o
	$(FC) $(FFLAGS) -c $< -o $@

src/matrices.o: src/matrices.f90 src/types.o src/potential.o src/internal_coordinates.o
	$(FC) $(FFLAGS) -c $< -o $@

src/eigenproblem.o: src/eigenproblem.f90
	$(FC) $(FFLAGS) -c $< -o $@

src/modify_structure.o: src/modify_structure.f90 src/types.o src/internal_coordinates.o
	$(FC) $(FFLAGS) -c $< -o $@

src/metrics.o: src/metrics.f90 src/types.o
	$(FC) $(FFLAGS) -c $< -o $@

src/main.o: src/main.f90 src/types.o src/io.o src/internal_coordinates.o src/matrices.o \
            src/eigenproblem.o src/modify_structure.o src/metrics.o
	$(FC) $(FFLAGS) -c $< -o $@

test_alignment.o: test_alignment.f90 src/types.o src/io.o src/metrics.o
	$(FC) $(FFLAGS) -c $< -o $@

$(TEST_EXEC): $(TEST_OBJ) src/types.o src/io.o src/metrics.o
	$(FC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f src/*.o src/*.mod *.o *.mod $(EXEC) $(TEST_EXEC)

test: $(EXEC)
	cp tests/*.pdb .
	./$(EXEC)

test_alignment: $(TEST_EXEC)
	cp tests/*.pdb .
	./$(TEST_EXEC)

.PHONY: all clean test test_alignment