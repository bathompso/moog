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
		print >>of, "atmosphere  1\nmolecules   2\nlines       1\nflux/int    0\ndamping     2\nplot        90"
		of.close()

		jnk = subprocess.check_output("expect %s/automoog automoog%d.mod" % (library_path, n), shell=True)

		# Open output file
		df = open('automoog%d.out' % n, 'r')
		lines = df.read().splitlines()
		df.close()

		# Create blank lists for statistics
		feh, std, corr = [], [], []
		fe, tcorr, ncorr = 0, 0.0, 0

		# Loop through lines and save important values
		for l in range(3,len(lines)-1):
			# Skip lines we don't care about
			if lines[l].find('abundance') < 0 and lines[l].find('E.P.') < 0: continue
	
			# Initial element line
			if lines[l].find('abundance') > 20:
				if lines[l].split()[4] == 'Fe' and lines[l].split()[5] == 'I': fe = 2
				elif lines[l].split()[4] == 'Fe': fe = 1
				else: fe = 0
	
			# Abundance summary line
			if lines[l].find('abundance') < 10 and lines[l].find('abundance') >= 0:
				if fe > 0:
					feh.append(float(lines[l].split()[3]))
				if fe == 2:
					std.append(numpy.abs(float(lines[l].split()[7])))
		
			# Correlation slopes
			if lines[l].find('E.P.') >= 0 and lines[l].find('correlation') >= 0 and lines[l+1].find('correlation') >= 0:
				if fe == 2:
					ncorr = 2
					tcorr = numpy.abs(float(lines[l].split()[11])) + numpy.abs(float(lines[l+1].split()[11]))
					if len(lines[l-1].split()) >= 12:
						ncorr += 1
						tcorr += numpy.abs(float(lines[l-1].split()[11]))
					corr.append(tcorr/ncorr)
		
		if len(corr) > 0: mcorr = numpy.mean(corr)
		else: mcorr = 0
		
		if len(std) > 0: ms = numpy.mean(std)
		else: ms = 0
		
		results.append([max(feh) - min(feh), ms, mcorr])
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

Dfeh, avg_std, corr_coeff = np.zeros(len(modfiles)), np.zeros(len(modfiles)), np.zeros(len(modfiles))
for r in range(len(results)):
	res = results[r]()
	for x in range(len(res)):
		Dfeh[x*len(results)+r] = res[x][0]
		avg_std[x*len(results)+r] = res[x][1]
		corr_coeff[x*len(results)+r] = res[x][2]
	
# Determine residual "score"
good = [x for x in range(len(modfiles)) if avg_std[x] < 0.05 and Dfeh[x] < 0.15]
print("%4s %4s %4s %4s      %6s  %6s  %6s" % ('Teff', 'logg', 'feh', 'vmic', 'DFe/H', 'std.', 'Corr.'))
for m in good:
	# Determine model parameters
	teff = float(modfiles[m][1:5])
	logg = float(modfiles[m][6:8])/10.0
	modfeh = float(modfiles[m][9:12])/100.0
	vmicro = float(modfiles[m][13:15])/10.0	
	print("%4d %4.1f %4.2f %4.1f      %6.3f  %6.3f  %6.3f" % (teff, logg, modfeh, vmicro, Dfeh[m], avg_std[m], corr_coeff[m]))
	
of = open("automoog.out", 'w')
print("%4s %4s %4s %4s      %6s  %6s  %6s" % ('Teff', 'logg', 'feh', 'vmic', 'DFe/H', 'std.', 'Corr.'), file=of)
order = np.argsort(avg_std)
for m in order:
	# Determine model parameters
	teff = float(modfiles[m][1:5])
	logg = float(modfiles[m][6:8])/10.0
	modfeh = float(modfiles[m][9:12])/100.0
	vmicro = float(modfiles[m][13:15])/10.0	
	print("%4d %4.1f %4.2f %4.1f      %6.3f  %6.3f  %6.3f" % (teff, logg, modfeh, vmicro, Dfeh[m], avg_std[m], corr_coeff[m]), file=of)
of.close()
	
	