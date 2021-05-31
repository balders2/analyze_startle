



import glob, os, sys, subprocess, shutil
for fn in glob.glob("*.acq"):
	base=fn.rsplit( ".", 1 )[ 0 ]
	matfn = base + '.mat'
	cmdn = 'acq2mat' + ' ' + fn + ' ' + matfn
	os.system(cmdn)
	subprocess.call(['./analyze_startle.exe', matfn])
	shutil.move("./" + fn, "./startle/" + fn)

	
for edffn in glob.glob("*.EDF"):
	subprocess.call(['./analyze_startle.exe', edffn])
	shutil.move("./" + edffn, "./startle/" + edffn)
