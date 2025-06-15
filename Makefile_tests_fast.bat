type src\matrix.cpp^
 src\AccelPointMass.cpp^
 src\Cheb3D.cpp^
 src\EccAnom.cpp^
 src\Frac.cpp^
 src\MeanObliquity.cpp^
 src\Mjday.cpp^
 src\Mjday_TDB.cpp^
 src\Position.cpp^
 src\R_x.cpp^
 src\R_y.cpp^
 src\R_z.cpp^
 src\sign_.cpp^
 src\timediff.cpp^
 src\AzElPa.cpp^
 src\IERS.cpp^
 src\Legendre.cpp^
 src\NutAngles.cpp^
 src\TimeUpdate.cpp^
 src\global.cpp^
 src\AccelHarmonic.cpp^
 src\EqnEquinox.cpp^
 src\JPL_Eph_DE430.cpp^
 src\LTC.cpp^
 src\NutMatrix.cpp^
 src\PoleMatrix.cpp^
 src\PrecMatrix.cpp^
 src\gmst.cpp^
 src\gast.cpp^
 src\MeasUpdate.cpp^
 src\G_AccelHarmonic.cpp^
 src\GHAMatrix.cpp^
 src\Accel.cpp^
 src\VarEqn.cpp^
 src\DEInteg.cpp^
 tests\tests.cpp^
 > tests\result.cpp  2>nul

g++ tests/result.cpp -lm -std=c++23 -o bin/tests.exe
del tests\result.cpp
cd bin
tests.exe
pause