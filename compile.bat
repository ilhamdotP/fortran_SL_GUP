gfortran -c tov-1.3.f
gfortran -o run.exe tov-1.3.o
del tov-1.3.o
@pause
run.exe
@pause
del run.exe