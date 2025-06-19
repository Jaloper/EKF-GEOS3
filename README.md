# Proyecto TTI 
## Instrucción de compilación aplicación principal:
g++ tests/EKF_GEOS3.cpp src/*cpp -lm -std=c++23 -o bin/main.exe
## Instrucción de compilación de tests unitarios:
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe