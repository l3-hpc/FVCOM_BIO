# BIO_TP
The TP bio model for FVCOM.

## More than one TP
(I am still testing this...)
To run more than one TP scenario, e.g. different sinking rates, and whether it sinks out:

TP.in
```
4           !Number of vars
"TP" "mg/L" !Name and units
0.00        !Minimum value for TP
6.0e-7      !Sinking rate, m/s
0           !1=sink out, 0=no sink out
"TP1" "mg/L" !Name and units
0.00        !Minimum value for TP
6.0e-7      !Sinking rate, m/s
1           !1=sink out, 0=no sink out
"TP2" "mg/L" !Name and units
0.00        !Minimum value for TP
6.0e-6      !Sinking rate, m/s
0           !1=sink out, 0=no sink out
"TP3" "mg/L" !Name and units
0.00        !Minimum value for TP
6.0e-6      !Sinking rate, m/s
1           !1=sink out, 0=no sink out
```

add more TPs to river and restart (use 'append'):
```
module load cpu/0.15.4 gcc/10.2.0 openmpi/4.0.4 nco/4.9.3 netcdf-c/4.7.4

ncks -O 2013_leem_fine_river_data.nc save_river_data.nc
ncap2 -s"TP1=TP" -v 2013_leem_fine_river_data.nc 2013_leem_fine_river_data.nc 
ncap2 -s"TP2=TP" -v 2013_leem_fine_river_data.nc 2013_leem_fine_river_data.nc 
ncap2 -s"TP3=TP" -v 2013_leem_fine_river_data.nc 2013_leem_fine_river_data.nc

ncks -O restart_tp.nc save_restart_tp.nc
ncap2 -s"TP1=TP" -v restart_tp.nc restart_tp.nc
ncap2 -s"TP2=TP" -v restart_tp.nc restart_tp.nc
ncap2 -s"TP3=TP" -v restart_tp.nc restart_tp.nc
```

Check the metadata:
```
ncdump -h 2013_leem_fine_river_data.nc 
ncdump -h restart_tp.nc 
```

Also do:
```
cp NUTRIENT_INI_1.dat NUTRIENT_INI_2.dat 
cp NUTRIENT_INI_1.dat NUTRIENT_INI_3.dat 
cp NUTRIENT_INI_1.dat NUTRIENT_INI_4.dat 
```

