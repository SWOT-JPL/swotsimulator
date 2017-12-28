import glob, os, sys, shutil
test_files = glob.glob('*.py')
test_files.remove('run_all.py')
test_files.remove('myspectools.py')
test_files.remove('compute_spectrum_geophys.py')
test_files.remove('plot_nadir.py')
test_files.remove('params.py')
'''
if py_path is None:
    py_path = '.'
else:
    py_path = os.pathsep.join(['.',py_path])
os.environ['PYTHONPATH'] = py_path
'''
if  (len(sys.argv) < 2):
    file_param='../example'+os.sep+'params_example.py'
    print("no params file specified, default is " +file_param)
else:
    file_param=str(sys.argv[1])
if os.path.isfile(file_param):
        shutil.copyfile(file_param, 'params.py')
else:
        print("Error: No such file: '%s'" % file_param)
        sys.exit()

for f in test_files:
    sys.stdout.write( "**********************************************\n")
    ff = os.path.join(sys.path[0],f)
    args = [sys.executable,ff]
    sys.stdout.write("Running %s\n" % f)
    status = os.spawnve(os.P_WAIT,sys.executable,args,os.environ)
    if status:
        sys.stdout.write('TEST FAILURE (status=%s)\n' % (status))

