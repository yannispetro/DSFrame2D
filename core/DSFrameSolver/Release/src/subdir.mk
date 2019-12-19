################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/1_read_input.f90 \
../src/2_Stiffness_Matrix.f90 \
../src/3_modifications.f90 \
../src/41_eigen.f90 \
../src/4_solver.f90 \
../src/5_element_forces.f90 \
../src/6_write_output.f90 \
../src/_MAIN.f90 \
../src/datmod.f90 \
../src/matrixtools.f90 

OBJS += \
./src/1_read_input.o \
./src/2_Stiffness_Matrix.o \
./src/3_modifications.o \
./src/41_eigen.o \
./src/4_solver.o \
./src/5_element_forces.o \
./src/6_write_output.o \
./src/_MAIN.o \
./src/datmod.o \
./src/matrixtools.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O3 -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/1_read_input.o: ../src/1_read_input.f90 src/datmod.o

src/2_Stiffness_Matrix.o: ../src/2_Stiffness_Matrix.f90 src/datmod.o

src/3_modifications.o: ../src/3_modifications.f90 src/datmod.o

src/41_eigen.o: ../src/41_eigen.f90 src/datmod.o

src/4_solver.o: ../src/4_solver.f90 src/datmod.o

src/5_element_forces.o: ../src/5_element_forces.f90 src/datmod.o

src/6_write_output.o: ../src/6_write_output.f90 src/datmod.o

src/_MAIN.o: ../src/_MAIN.f90 src/datmod.o

src/datmod.o: ../src/datmod.f90

src/matrixtools.o: ../src/matrixtools.f90


