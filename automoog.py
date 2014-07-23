from __future__ import print_function, division
import subprocess, sys, pp
import numpy as np
from time import time

##### CHANGE VARIABLES HERE #####
library_path = '/automoog_library'
#################################

def modelfit(modfiles, ncpus, n, linefile, library_path):
	results = []
	for m in range(len(modfiles)):
		if m % ncpus != n: continue

		# Create new model file
		of = open("automoog%d.mod" % n, 'w')
		print >>of, "abfind"
		print >>of, "terminal  'x11'"
		print >>of, "standard_out  'automoog%d'" % (n)
		print >>of, "summary_out  'automoog%d.out'" % (n)
		print >>of, "model_in  '%s/%s'" % (library_path, modfiles[m])
		print >>of, "lines_in  '%s'" % (linefile)
		print >>of, "atmosphere  1\nmolecules   2\nlines       1\nflux/int    0\ndamping     2\nplot        99"
		of.close()

		jnk = subprocess.check_output("expect %s/automoog.scr automoog%d.mod" % (library_path, n), shell=True)

		# Open output file
		df = open('automoog%d.out' % n, 'r')
		lines = df.read().splitlines()
		df.close()

		# Create blank lists for statistics
		feh, std, b = [], [], []
		fe = 0

		# Loop through lines and save important values
		for l in range(3,len(lines)-1):
			# Skip lines we don't care about
			if lines[l].find('abundance') < 0 and lines[l].find('E.P.') < 0: continue
	
			# Initial element line
			if lines[l].find('abundance') > 20:
				if lines[l].split()[4] == 'Fe': fe = 1
				else: fe = 0
	
			# Abundance summary line
			if lines[l].find('abundance') < 10 and lines[l].find('abundance') >= 0:
				if fe == 1:
					feh.append(float(lines[l].split()[3]))
					std.append(float(lines[l].split()[7]))
		
			# Correlation slopes
			if lines[l].find('E.P.') >= 0 and lines[l+1].find('correlation') >= 0:
				b.append(float(lines[l].split()[7]) - float(lines[l+1].split()[7]))
		
		results.append([max(feh) - min(feh), numpy.mean(std), numpy.mean(b)])
	return results

linefile = sys.argv[1]

# Find all usable models
modfiles = subprocess.check_output("ls %s | grep v" % (library_path), shell=True).splitlines()
	
# Create job server
jobServer = pp.Server()

# Run everything
start_time = time()
results = [jobServer.submit(modelfit, (modfiles, jobServer.get_ncpus(), n, linefile, library_path), (), ("numpy","subprocess")) for n in range(jobServer.get_ncpus())]
for r in results: tmp = r()
end_time = time()
print("ELAPSED: %d min." % ((end_time - start_time)/60))
subprocess.check_output("ls automoog* | grep -v py | xargs rm", shell=True)

Dfeh, avg_std, Db = np.zeros(len(modfiles)), np.zeros(len(modfiles)), np.zeros(len(modfiles))
for r in range(len(results)):
	res = results[r]()
	for x in range(len(res)):
		Dfeh[x*len(results)+r] = res[x][0]
		avg_std[x*len(results)+r] = res[x][1]
		Db[x*len(results)+r] = res[x][2]
	
# Determine residual "score"
score = np.abs(Dfeh) + np.abs(avg_std) + np.abs(Db)
sort_score = np.argsort(score)
print("%4s %4s %4s %4s      %6s  %6s  %6s" % ('Teff', 'logg', 'feh', 'vmic', 'Fe/H', 'Sigma', 'Intcp'))
for m in sort_score[0:20]:
	# Determine model parameters
	teff = float(modfiles[m][1:5])
	logg = float(modfiles[m][6:8])/10.0
	modfeh = float(modfiles[m][9:12])/100.0
	vmicro = float(modfiles[m][13:15])/10.0	
	print("%4d %4.1f %4.2f %4.1f      %6.3f  %6.3f  %6.3f" % (teff, logg, modfeh, vmicro, Dfeh[m], avg_std[m], Db[m]))
	
	