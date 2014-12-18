rm -f get_diag.so
f2py -c -m get_diag diag.f90 -L/usr/lib -lblas -llapack
