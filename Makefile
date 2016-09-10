# Makefile : ITPRPRP
#
# TESTED
#
#    GNU Make 3.81
#    GNU Fortran (Homebrew gcc 5.3.0) 5.3.0
#
# LAST TESTED
#
#    GNU Make 3.81
#    GNU Fortran (Homebrew gcc 5.3.0) 5.3.0
#    2016/09/09

COMPILER          := gfortran 
COMPILER_OPTIONS  := -fimplicit-none -ffree-form -ffree-line-length-none -std=gnu -O2 -mtune=native
SOURCE_DIRECTORY  := source
LAPACK_DIRECTORY  := libraries/lapack
BUILD_DIRECTORY   := build
TARGET_EXECUTABLE := itpprp.x 
SOURCE            := $(shell find $(SOURCE_DIRECTORY) -type f -name *.f)
OBJECTS           := $(patsubst $(SOURCE_DIRECTORY)/%,$(BUILD_DIRECTORY)/%,$(SOURCE:.f=.o)) 

$(TARGET_EXECUTABLE): $(OBJECTS) $(BUILD_DIRECTORY)/zgtsv.o $(BUILD_DIRECTORY)/xerbla.o
	@echo " Linking..."
	@echo " $(COMPILER) $^ -o $(TARGET_EXECUTABLE)"; $(COMPILER) $^ -o $(TARGET_EXECUTABLE)

$(BUILD_DIRECTORY)/%.o: $(SOURCE_DIRECTORY)/%.f
	@mkdir -p $(BUILD_DIRECTORY)
	@echo " $(COMPILER) $(COMPILER_OPTIONS) -c -o $@ $<"; $(COMPILER) $(COMPILER_OPTIONS) -c -o $@ $<

$(BUILD_DIRECTORY)/zgtsv.o: $(LAPACK_DIRECTORY)/zgtsv.f
	$(COMPILER) $(COMPILER_OPTIONS) -o $(BUILD_DIRECTORY)/zgtsv.o -c $(LAPACK_DIRECTORY)/zgtsv.f 

$(BUILD_DIRECTORY)/xerbla.o: $(LAPACK_DIRECTORY)/xerbla.f
	$(COMPILER) $(COMPILER_OPTIONS) -o $(BUILD_DIRECTORY)/xerbla.o -c $(LAPACK_DIRECTORY)/xerbla.f

.PHONY: clean
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILD_DIRECTORY) $(TARGET_EXECUTABLE)"; $(RM) -r $(BUILD_DIRECTORY) $(TARGET_EXECUTABLE)
