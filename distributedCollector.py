# CartRing Distributed Result Collector
# Andrew Stershic - DCML
# ajs84@duke.edu
# 5/31/2013

# Import necessary packages
#from subprocess import call
import time
import os
import sys
import string

# --------------------------------
# Main function
#print "LAUNCHING DISTRIBUTED COLLECTOR"
argStr = str(sys.argv)
argStr = argStr.translate(string.maketrans('', ''), '!@#$[]\',')
argList = argStr.split()

if len(argList) >= 4:
	base = argList[1]
	dest = argList[2]
	hostList = argList[3:]
else:
	base = ""
	dest = ""
	hostList=[]
	print "***Empty hostlist"
	print argStr
	
r = dest.rfind("results")
dest = dest[0:r-1]

command = dest + "/./clean.sh"
os.system(command)		

for host in hostList:
	command1 = "scp -r " + host + ":" + base + " " + dest 
	command2 = "ssh " + host + " " + base + "/./clean.sh"
	command3 = "ssh " + host + " rm " + base + "/*.log"
	print command1
	os.system(command1)
	#print command2
	os.system(command2)
	os.system(command3)
	
	
print "*** results collected from compute nodes"