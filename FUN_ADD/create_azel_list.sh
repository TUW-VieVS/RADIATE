# This batch-file creates an executable file for the creation of an azel file list:
gfortran -c ../FUN_TEXT/module_date_type_definition.f90
gfortran -c ../FUN_TEXT/module_date2mjd.f90
gfortran -c ../FUN_TEXT/module_mjd2date.f90 
gfortran create_azel_name_list.f90 module_date2mjd.o module_mjd2date.o module_date_type_definition.o -o azel_list

 


