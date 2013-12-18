# CartRing Automated Runs
# Andrew Stershic - DCML
# ajs84@duke.edu
# 4/17/2013

# Import necessary packages
from subprocess import call
import time
import os

# Define tail function which will return the last 'window' lines
def tail( fh, window=1 ):
	f = open(fh, 'r')
	try: 
		BUFSIZ = 1024
		f.seek(0, 2)
		bytes = f.tell()
		size = window
		block = -1
		data = []
		while size > 0 and bytes > 0:
		    if (bytes - BUFSIZ > 0):
		        # Seek back one whole BUFSIZ
		        f.seek(block*BUFSIZ, 2)
		        # read BUFFER
		        data.append(f.read(BUFSIZ))
		    else:
		        # file too small, start from begining
		        f.seek(0,0)
		        # only read what was not read
		        data.append(f.read(bytes))
		    linesFound = data[-1].count('\n')
		    size -= linesFound
		    bytes -= BUFSIZ
		    block -= 1
		return '\n'.join(''.join(data).splitlines()[-window:])
	finally:
		f.close()
		
def clearMultiHost(hostList,base):
	for host in hostList:
		command2 = "ssh " + host + " " + base + "/./clean.sh"
		command3 = "ssh " + host + " rm " + base + "/*.log"
		os.system(command2)
		os.system(command3)


# --------------------------------
# Main function

# Get to the right place
#base = "/Users/andrewstershic/Code/axisymmetricRing/"
#base = "/home/ajs84/Code/axisymmetricRing/"
base = os.getcwd()
call(["clear"])
os.chdir(base)

# Clean out old results (if needed)
if not os.path.exists("results/automatedRuns"):
    os.makedirs("results/automatedRuns")
filelist = [ f for f in os.listdir("results/automatedRuns") if f.endswith(".dat") ]
for f in filelist:
    os.remove("results/automatedRuns/" + f)

# Compile Program
#call(["make","clean"])
#call(["make"])

# Initialize run series information
# run <n> independent problems per node
n = 1
# label the run series, default ""
series = ""
# spread the load over <p> processors per node, max = 8
p = 8
# run on <m> nodes, max = 16. Total number of processors = <m*p>
m = 10
hostList = ["f1","f2","f3","f4","f5","f6","f7","f8","f9","f10"]
assert(len(hostList) == m)

print ""
print "---------------------------------------------"
print "**************BEGIN EXECUTION****************"
print "---------------------------------------------"
print ""
	
# Start the timer	
#start = time.clock()
startT = time.time()

# Run the series
for x in range(n):
	print ""
	print "Execution number %d" % x
	if (p == 1):
		call(["./bin/main.exe"])
	else:
		
		clearMultiHost(hostList,"/tmp/results")
		
		hostStr = ""
		hostStr2 = ""
		for host in hostList:
			hostStr += host + " "
			hostStr2 += host + ","
		hostStr2 = hostStr2[0:-1]
		os.system("mpirun -H " + hostStr2 + " -np " + str(m*p) + " ./bin/main.exe " +  hostStr)
		#call(["mpirun","-np",str(p),"./bin/main.exe " + hostList])
	filelist = [ f for f in os.listdir("results/datFiles") if f.endswith(".dat") ]
	for f in filelist:
		f_raw = f.replace(".dat","")
		fileName = "results/automatedRuns/" + series + f_raw + "-" + str(x)  + ".dat"
		print fileName
		os.rename("results/datFiles/"+f,fileName)

# Show timing results
#elapsed = (time.clock() - start)
elapsedT = (time.time() - startT)
print ""
#print "Average run-time = ",elapsed/n, " seconds"
print "Average run-time = ",elapsedT/n, " seconds"

print ""
print "---------------------------------------------"
print "***************END EXECUTION*****************"
print "---------------------------------------------"
print ""

# Data Collection and averaging

# Need to collect and average data
min_frag_size = 0
mean_frag_size = 0;
num_frags = 0
for x in range(n):
	f = "results/automatedRuns/" + series + "fraginfo-" + str(x) + ".dat"
#	f = os.popen("results/automatedRuns/fraginfo-" + str(x) + ".dat")
	line = tail(f,1)
	line = line.split()
	min_frag_size += float(line[6])/n
	mean_frag_size += float(line[3])/n
	num_frags += float(line[1])/n

fragLengths = []
for x in range(n):
	f = "results/automatedRuns/" + series + "fraghisto-" + str(x) + ".dat"
	fh = open(f)
	lines = fh.readlines()
	fh.close()

	foundSizes = False
	while ((foundSizes == False) and (len(lines) > 0)):
		if "Sizes:" in lines[0]:
			foundSizes = True
			lines.pop(0)
		else:
			lines.pop(0)

	for y in lines:	
		fragLengths.append(y.rstrip())

fragfile = open("results/automatedRuns/" + series + "allFragmentLengths.dat", 'w')
for item in fragLengths:
	print>>fragfile, item
fragfile.close()


print 'Average minimum fragment size = ', min_frag_size
print 'Average mean fragment size = ', mean_frag_size
print 'Average number of fragments = ', num_frags
