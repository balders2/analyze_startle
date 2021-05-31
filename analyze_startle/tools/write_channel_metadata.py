



import glob, os, sys, subprocess, shutil
for fn in glob.glob("*.acq"):
	base=fn.rsplit( ".", 1 )[ 0 ]
	matfn = base + '.mat'
	cmdn = 'acq2mat' + ' ' + fn + ' ' + matfn
	os.system(cmdn)
	subprocess.call(['./write_channel_metadata.exe', matfn])
	shutil.move("./" + fn, "./startle/" + fn)