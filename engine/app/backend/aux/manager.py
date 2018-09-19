import os
import sys
import json
import getopt

#
# External variables
#

opts, _ = getopt.getopt(sys.argv[1:],'s:l:',['suffix=','label='])

for opt, arg in opts:
	if opt in ('-l','--label'):
		label = arg
	if opt in ('-s','--suffix'):
		suffix = arg

#
# Print specific label
#

with open('../../static/results/config-'+suffix+'.json','r') as fj:
	content = json.load(fj)

print(content[label])