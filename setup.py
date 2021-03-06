import os
import sys
import fileinput
import platform
import subprocess

from platform import system as current_platform
from shutil import copy, move, copytree
from glob import glob

__version__ = '2013.02'

# We need to build MOOG and MOOGSILENT before they get moved to the scripts/
# directory so that they can be moved into the $PATH
if 'build' in sys.argv:

	# Identify the platform
	platform = current_platform()

	# Check for platform first
	if platform not in ('Darwin', 'Linux'):
		sys.stderr.write("Platform '%s' not recognised!\n" % platform)
		sys.exit()


	if os.getuid() != 0:
		sys.stderr.write("Permission denied: Sudo access is required!\n")
		sys.exit()

	# We are sudo; with great power comes great responsibility.

	# By default, we will use 32bit 
	is_64bits = False

	# Which system are we on?
	if platform == 'Darwin':
		run_make_files = ('Makefile.mac', 'Makefile.macsilent')
		machine = 'mac'

	elif platform == 'Linux':

		machine = 'pcl'
		is_64bits = sys.maxsize > 2**32

		if is_64bits:
			run_make_files = ('Makefile.rh64', 'Makefile.rh64silent')

		else:
			run_make_files = ('Makefile.rh', 'Makefile.rhsilent')


	# Check for gfortran or g77
	def system_call(command):
		""" Perform a system call with a subprocess pipe """
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
		return process.communicate()[0]

	# Look for g77 and gfortran
	g77_exists = len(system_call("which g77")) > 0
	gfortran_exists = len(system_call("which gfortran")) > 0

	# If we have the choice, use gfortran
	if gfortran_exists:

		if is_64bits:
			fortran_vars = "FC = gfortran -m64\nFFLAGS = -Wall -O4 -ffixed-line-length-72 -ff2c"

		else:
			fortran_vars = "FC = gfortran\nFFLAGS = -Wall -O4 -ffixed-line-length-72 -ff2c"

	elif g77_exists:
		
		if platform == 'Linux':
			fortran_vars = 'FC = g77 -Wall'

		else:
			fortran_vars = 'FC = g77 -w'
		
	else:
		sys.stderr.write("Could not find g77 or gfortran on the system!\n")
		sys.exit()

	# Get our directories relative to the current path
	repository_dir = os.path.dirname(os.path.realpath(__file__))

	# We need a moog data directory
	data_dir = os.path.expanduser('/.moog')
	if not os.path.exists(data_dir):
		system_call('mkdir %s' % data_dir)

	# Copy files to data directory
	src_dir = os.path.join(repository_dir, 'src')
	data_files = glob('%s/*.dat' % src_dir)
	[copy(data_file, '%s/%s' % (data_dir, os.path.basename(data_file), )) for data_file in data_files]
	
	aqlib = "AQLIB = %s" % os.path.join(repository_dir, 'lib/aqlib')
	smlib = "SMLIB = %s" % os.path.join(repository_dir, 'lib/smlib')

	configuration = "\n".join([fortran_vars, aqlib, smlib])

	# Update the makefiles with the proper SMLIB and AQLIB
	run_make_files = [os.path.join(repository_dir, 'src', filename) for filename in run_make_files]
	hardcoded_moog_files = [os.path.join(repository_dir, 'src', filename) for filename in ('Begin.f', 'Moog.f', 'Moogsilent.f')]

	# Setup: Move and create copies of the original
	for make_file in run_make_files:
		move(make_file, make_file + '.original')
		copy(make_file + '.original', make_file)

	for moog_file in hardcoded_moog_files:
		move(moog_file, moog_file + '.original')
		copy(moog_file + '.original', moog_file)

	# Update the run make files with the configuration
	for line in fileinput.input(run_make_files, inplace=True):
		line = line.replace('#$CONFIGURATION', configuration)

		sys.stdout.write(line)

	# Update the MOOG files
	for line in fileinput.input(hardcoded_moog_files, inplace=True):
		line = line.replace('$SRCDIR', src_dir)
		line = line.replace('$DATADIR', data_dir)
		line = line.replace('$MACHINE', machine)

		sys.stdout.write(line)

	# Run the appropriate make files
	for make_file in run_make_files:
		os.system('cd src;make -f %s' % make_file)

	# Cleanup files: Replace with original files
	[move(moog_file + '.original', moog_file) for moog_file in hardcoded_moog_files if os.path.exists(moog_file + '.original')]
	[move(make_file + '.original', make_file) for make_file in run_make_files if os.path.exists(make_file + '.original')]

	# Copy the AquaTerm framework
	if not os.path.exists('/Library/Frameworks/AquaTerm.framework/'):
		try:
			system_call('cp -R %s /Library/Frameworks/AquaTerm.framework/' % os.path.join(repository_dir, 'lib/AquaTerm.framework/'))

		except:
			sys.stdout.write("AquaTerm framework could not be installed to /Library/Frameworks/AquaTerm.framework\n")

		else:
			sys.stdout.write("AquaTerm framework copied to /Library/Frameworks/AquaTerm.framework\n")


# If binaries are already built, we just need to install them to the right place
elif 'install' in sys.argv:
	# Define all the necessary paths
	repository_dir = os.path.dirname(os.path.realpath(__file__))
	data_dir = os.path.expanduser('/.moog')

	# Make the data directory to hold data files
	if not os.path.exists(data_dir):
		system_call('mkdir %s' % data_dir)
		
	# Move the data files to the data directory
	os.system('cp %s/src/Barklem* %s/' % (repository_dir, data_dir))
	
	# Move the pre-compiled binaries to the right location
	os.system('cp %s/bin/* /usr/local/bin/' % repository_dir)
	
	# Move the necessary frameworks to the right location
	if not os.path.exists('/usr/local/lib/libg2c.0.dylib'):
		os.system('cp %s/lib/libg2c.0.dylib /usr/local/lib/' % repository_dir)
	
	# Copy the AquaTerm framework
	if not os.path.exists('/Library/Frameworks/AquaTerm.framework/'):
		try: system_call('cp -R %s /Library/Frameworks/AquaTerm.framework/' % os.path.join(repository_dir, 'lib/AquaTerm.framework/'))
		except: sys.stdout.write("AquaTerm framework could not be installed to /Library/Frameworks/AquaTerm.framework\n")
		else: sys.stdout.write("AquaTerm framework copied to /Library/Frameworks/AquaTerm.framework\n")
