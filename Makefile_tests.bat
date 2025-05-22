g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe
cd bin
tests.exe
pause