from __future__ import print_function, division
import numpy as np
import subprocess, sys, os

# Specify output directory from command line
outdir = sys.argv[1]

# Check to see whether output directory exists, or can be created
if not os.path.exists(outdir):
	try:
		subprocess.check_output("mkdir %s" % outdir, shell=True)
	except:
		print("Could not create directory '%s'. Aborting..." % outdir)
		exit(0)

# Move automoog expect script to correct location
subprocess.check_output("cp automoog.scr %s/" % (outdir), shell=True)

# Compile the makekurucz executable
subprocess.check_output("gfortran src/makekurucz3.f src/phdlibsd.f -ffixed-line-length-132 -o makekurucz", shell=True)

teffs = np.arange(5200, 6510, 10)
loggs = np.arange(3.0, 5.6, 0.1)
fehs = [0]
vmicros = [1.8]

for t in teffs:
	for g in loggs:
		for f in fehs:
			for v in vmicros:
				of = open("tmp.in", 'w')
				print(t, g, f, v, file=of)
				print("NOVER", file=of)
				print("", file=of)
				of.close()
				
				try:
					subprocess.check_output("./makekurucz < tmp.in", shell=True)
					jnk = subprocess.check_output("mv FINALMODEL %s/t%4dg%02df%03dv%02d" % (outdir, t, g*10, f*100, v*10), shell=True)
				except:
					pass

subprocess.check_output("rm tmp.in; rm MOD*", shell=True)
