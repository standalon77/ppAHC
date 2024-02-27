################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/ClientSocket.cpp \
../src/PPkNNCSPMain.cpp \
../src/PaillierCrypto.cpp \
../src/Socket.cpp 

CPP_DEPS += \
./src/ClientSocket.d \
./src/PPkNNCSPMain.d \
./src/PaillierCrypto.d \
./src/Socket.d 

OBJS += \
./src/ClientSocket.o \
./src/PPkNNCSPMain.o \
./src/PaillierCrypto.o \
./src/Socket.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/ClientSocket.d ./src/ClientSocket.o ./src/PPkNNCSPMain.d ./src/PPkNNCSPMain.o ./src/PaillierCrypto.d ./src/PaillierCrypto.o ./src/Socket.d ./src/Socket.o

.PHONY: clean-src

