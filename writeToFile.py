import numpy as np
import uproot

# ----------------------------------------------------------------------------------------------
# Open data file
fIn = uproot.open('S2-tuple_coins.root')

# Select trees of interest
clusterTree = fIn['DefaultReconstruction/clusters']
muonTree = fIn['DefaultReconstruction/muons']
nsEsCoinsTree = fIn['DefaultReconstruction/nsEsClusCoins']

# Select leaves of interest
clusterTreeEntries=clusterTree.arrays(['start', 'rejected', 'MuonTag', 'NSTag', 'ESTag', 'MuonTreePos', 'chanID', 'chanActive', 'cubes', 'wfTreePos', 'chanAmplitude', 'onionLayer'])
muonTreeEntries = muonTree.arrays(['htime','startCube','endCube','type','stopped'])
nsEsCoinsEntries = nsEsCoinsTree.arrays(['time','delt','delx','dely','delz','delr','promptEnergy','NSClusTreePos','ESClusTreePos'])

# -----------------------------------------------------------------------------------------------
# Check for muon decay with michel electron:
def michelCheck(tree):
	event1, event2 = [], []
	for iEntry in range(len(tree) - 1):
		if clusterTreeEntries['MuonTag'][iEntry]:
			# position in muonTree of event
			muonPos = clusterTreeEntries['MuonTreePos'][iEntry]

			# Determine time of muon events
			time_old = muonTreeEntries['htime'][muonPos] * 1e-3 # microseconds!
			time = muonTreeEntries['htime'][muonPos + 1] * 1e-3

			delta_t = abs(time - time_old)

			# determine endpoint of muon track and startpoint of the next muon track
			end = muonTreeEntries['endCube'][muonPos]
			start = muonTreeEntries['startCube'][muonPos + 1]

			""" 
			Note: The distance is measured in cubes, Since the cubes have dimensions (5cm x 5cm x 5cm), we can easily convert the units to cm. 
			"""
			delta_x = abs(start[0] - end[0]) * 5
			delta_y = abs(start[1] - end[1]) * 5
			delta_z = abs(start[2] - end[2]) * 5

			delta_r = np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)

			# Check if the particle stopped within the detector
			stopped = muonTreeEntries['stopped'][muonPos]

			if (delta_t <= 6) and (delta_r <= 30) and stopped:
				event1.append(iEntry)
				event2.append(iEntry + 1)

	return event1, event2
	
# ----------------------------------------------------------------------------------------------
# check for IBD event (simplified):
def IBDcheck(tree):
	nsPos, esPos = [], []

	# IBD cuts
	for iEntry in range(len(tree)):
		# IBD cuts
		if(not(100000>=nsEsCoinsEntries['delt'][iEntry]>0)): continue
		if(not(2>=nsEsCoinsEntries['delx'][iEntry]>=-2)): continue
		if(not(2>=nsEsCoinsEntries['dely'][iEntry]>=-2)): continue
		if(not(3>=nsEsCoinsEntries['delz'][iEntry]>=-2)): continue
		if(not(3>=nsEsCoinsEntries['delr'][iEntry]>0)): continue
		if(not(6.5>nsEsCoinsEntries['promptEnergy'][iEntry]>1)): continue
				
		nsPos.append(nsEsCoinsEntries['NSClusTreePos'][iEntry])
		esPos.append(nsEsCoinsEntries['ESClusTreePos'][iEntry])

	return nsPos, esPos
# ---------------------------------------------------------------------------------------------
# Write .txt file with selected events
def WriteToFile(criteria, filename, tree=clusterTreeEntries['start'], numberOfEventsTogether=2):
	# write .txt file with name filename
	file = open(filename, 'w')

	# write first line
	for i in range(numberOfEventsTogether):
		file.write('Event of Type {} \t'.format(i + 1))
	file.write('\n')

	# write subsequent lines
	if numberOfEventsTogether == 1:
		events = criteria(tree)
		for event in events:
			file.write('{}\n'.format(event))

	elif numberOfEventsTogether == 2:
		events1, events2 = criteria(tree)
		for i, event in enumerate(events1):
			file.write('{} \t {}\n'.format(event, events2[i]))
	file.close()

# ---------------------------------------------------------------------------------------------
# test the written function:
WriteToFile(criteria=michelCheck, filename='michelEvents.txt')
WriteToFile(criteria=IBDcheck, filename='IBDevents.txt', tree=nsEsCoinsEntries['time'], numberOfEventsTogether=2)
