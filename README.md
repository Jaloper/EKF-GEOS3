# EKF-GEOS3: C++ Extended Kalman Filter for GEOS-3 Orbit Determination

## Build Instructions

### Compile Main Application
```bash
g++ tests/EKF_GEOS3.cpp src/*cpp -lm -std=c++23 -o bin/main.exe
```
### Compile & Run Unit Tests
```bash
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe
```
