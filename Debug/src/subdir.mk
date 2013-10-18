################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BranchAndCut.cpp \
../src/ColumnGeneration.cpp \
../src/Dinic.cpp \
../src/Edge.cpp \
../src/Instance.cpp \
../src/Kruskals.cpp \
../src/LocalSearch.cpp \
../src/MinCut.cpp \
../src/PathRelinking.cpp \
../src/TreeGlover.cpp \
../src/Utils.cpp \
../src/VariableFixing.cpp \
../src/main.cpp 

OBJS += \
./src/BranchAndCut.o \
./src/ColumnGeneration.o \
./src/Dinic.o \
./src/Edge.o \
./src/Instance.o \
./src/Kruskals.o \
./src/LocalSearch.o \
./src/MinCut.o \
./src/PathRelinking.o \
./src/TreeGlover.o \
./src/Utils.o \
./src/VariableFixing.o \
./src/main.o 

CPP_DEPS += \
./src/BranchAndCut.d \
./src/ColumnGeneration.d \
./src/Dinic.d \
./src/Edge.d \
./src/Instance.d \
./src/Kruskals.d \
./src/LocalSearch.d \
./src/MinCut.d \
./src/PathRelinking.d \
./src/TreeGlover.d \
./src/Utils.d \
./src/VariableFixing.d \
./src/main.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/lapo/cplex/cplex125/concert/include -I/home/lapo/cplex/cplex125/cplex/include -O1 -g3 -Wall -c -fmessage-length=0 -m64  -fPIC -fexceptions -DIL_STD -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


