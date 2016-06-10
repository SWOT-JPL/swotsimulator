import sys, os
import swotsimulator.run_simulator as run_simulator

if  (len(sys.argv) < 2):
    file_param=os.getcwd()+os.sep+'example'+os.sep+'params_example_nadir.txt'

    print("no params file specified, default is " +file_param)
else:
    file_param=str(sys.argv[1])
run_simulator.run_nadir(file_param)
