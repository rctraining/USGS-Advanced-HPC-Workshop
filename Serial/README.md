# Serial optimization

## Building a program

### Pi directory

1. Build the program
```
make
make debug
```
2. Review optimization report
```
less Release/pi_serial.optrpt
```
3. Profiling with Allinea map
```
map Release/pi_serial
```
4. Compare profile with the debug version
5. Using vtune
```
bash /opt/intel/vtune-amplifier_xe/amplxe-vars.sh
amplxe-cl -collect Release/pi_serial < input.txt
amplxe-gui
```

### F90_matrix-multiply

Perform performance analysis of code

What is the major bottleneck?

### CPP_matrix-multiply

Same as fortran with index calculation in a 1D array

### CPP_boost_matrix_multiply

Same as Fortran 90 but with ```boost::multiarray```
