#!/bin/bash
make clean
make ib
echo "Output:==============================================================="
echo ""
##valgrind --tool=memcheck --leak-check=yes mpirun -np 3 ./bin/ib
mpirun -np 2 ./bin/ib -malloc_dump -ksp_monitor_true_residual -ksp_converged_reason  #-ksp_view -log_view #-on_error_attach_debugger gdb

