## Viterbi Decoder
# Nicholas McKibben
# ECEn 770
# 2018-04-14

from numpy import *

## State Object
class State(object):
	def __init__(self,thisState,nextStates,prevStates,output):
		self.this = thisState
		self.next = nextStates
		self.prev = prevStates
		self.output = output

class Branch(object):
	def __init__(self,states,metric=0):
		self.states = states
		self.metric = 0

class Occupied(object):
	def __init__(self):
		self.occupied_states = set()

	def get(self):
		return(self.occupied_states)

	def add(self,states):
		for s in states:
			self.occupied_states.add(s)

class Accepted(object):
	def __init__(self,t,prevState,thisState,metric):
		self.t = t
		self.cum_metric = metric
		self.this = thisState
		self.prev = prevState

def getKey(t,thisState,nextState):
	key = '%d %d %d' % (t,thisState,nextState)
	return(key)

def path_metric(a,b):
	return(sum(1 for x,y in zip(a,b) if x != y))

## Convolutional Encoder
def convencode(m,g):
    # Convolve for each g[j]
    c = zeros([ g.shape[0],(m.shape[0] + g.shape[1] -1) ])
    for ii in range(0,g.shape[0]):
        c[ii,:] = mod(convolve(m,g[ii,:]),2)
    
    # Interleave the rows
    c = reshape(c,[ 1,-1 ],order='F')
    return(c)

if __name__ == '__main__':

	## (1) Encoder
	# Write a computer program that performs the encoding operation
	# for a convolutional code with transfer matrix:
	g1 = array([ 1,0,1 ])
	g2 = array([ 1,1,1 ])

	# Test the encoder
	m = array([ 1,1,0,0,1,0,1 ])
	c = convencode(m,array([ g1,g2 ]))
	#print(c)

	## (2) Hard decision Viterbi
	# Program and test the hard-decision Viterbi decoder for this
	# encoder.  Your final result should be a BER performance plot
	# showing the performance of your code as a function of Eb/N0.
	# Include a curve that shows the theoretical performance for
	# uncoded BPSK. You may want to review Lab 1 material.  Do not
	# forget to incoprorate the rate of the code when generating
	# noise at specific SNR levels as you (hopefully) did in Lab1.
	# Make the y-axis log-scale.

	# Let's use Algorithm 12.1 in the book.
	## (1) Input
	r = array([ 1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,1 ])
	print('Start with r:')
	print(r.astype(float))

	## (3) Initialize
	# Set M(0) = and M(p) = Inf for p = 1,2,...,2^v - 1 to be
	# the initial path costs
	buff_len = 2**4
	M = array(inf*ones([ buff_len ]))
	M[0] = 0
	#print(M)

	# Get some states
	q1 = State(0,[ 0,1 ],[ 0,2 ],[ '00','11' ])
	q2 = State(1,[ 2,3 ],[ 0,2 ],[ '01','10' ])
	q3 = State(2,[ 0,1 ],[ 1,3 ],[ '11','00' ])
	q4 = State(3,[ 2,3 ],[ 1,3 ],[ '10','01' ])
	qs = [ q1,q2,q3,q4 ]
	#print(qs[1].output[1])

	# Set up occupied states
	currStates = Occupied()
	currStates.add([ 0 ]) # start in the 0th state

	# Set up accepted states array
	acceptedStates = dict()
	key = getKey(0,0,0)
	acceptedStates[key] = Accepted(0,0,0,0)

	# Set up the branches
	branches = []

	## (4) Set P = null for initial paths
	P = array(zeros([ buff_len ]))
	#print(P)

	## (5) Set t = 0
	t = 0
	## (6) Begin
	for idx in range(0,r.shape[0],2):
		## (7) For each state q at time t+1
		t = t + 1
		#print('Current r is %d%d, t = %d' % (r[idx],r[idx+1],t))
		curr_qs = currStates.get()
		for ii in list(curr_qs):
			## (8) Find the path metric for each path to state qi
			# There are two paths from each state corresponding
			# inputs 0,1
			currStates.add(qs[ii].next)
			for jj in [ 0,1 ]:
				metric = path_metric(''.join(str(e) for e in [ r[idx],r[idx+1] ]),qs[ii].output[jj])
				#print('%d -> %d is %d' % (qs[ii].this,qs[ii].next[jj],metric))

				# Add the path to the list of accepted states
				key = getKey(t,qs[ii].this,qs[ii].next[jj])
				prevKey = getKey(t-1,qs[ii].prev[0],qs[ii].this)
				if prevKey not in acceptedStates:
					# Try the other path
					prevKey = getKey(t-1,qs[ii].prev[1],qs[ii].this)

				if prevKey in acceptedStates:
					cum_metric = acceptedStates[prevKey].cum_metric + metric
				else:
					cum_metric = metric

				acceptedStates[key] = Accepted(t,qs[ii].this,qs[ii].next[jj],cum_metric)

		# Pick the branches with minimum metric
		#print('Summary of time step t = %d' % t)
		for cs in curr_qs:
			m = inf
			for jj in [ 0,1 ]:
				key = getKey(t,qs[cs].prev[jj],cs)
				if key in acceptedStates:
					if m >= acceptedStates[key].cum_metric:
						# Choose this metric!
						m = acceptedStates[key].cum_metric
						if jj is 1:
							del acceptedStates[getKey(t,qs[cs].prev[0],cs)]
		
			key = getKey(t,qs[cs].prev[0],cs)
			if key not in acceptedStates:
				key = getKey(t,qs[cs].prev[1],cs)
			#print('\t%d, metric = %d' % (cs,acceptedStates[key].cum_metric))

	# Trace back through the accepted states to find the recovered signal
	# Start at the last t we have
	#print('Traceback:')
	MLPath = []
	for ti in range(t,-1,-1):
		#print('\tt = %d' % ti)
		cost = inf
		MLkey = None
		for qi in qs:
			key = getKey(ti,qi.prev[0],qi.this)
			prev = qi.prev[0]
			if key not in acceptedStates:
				key = getKey(ti,qi.prev[1],qi.this)
				prev = qi.prev[1]
			if key in acceptedStates:
				#print('\t\t%d -> %d, metric = %d' % (prev,qi.this,acceptedStates[key].cum_metric))

				# find the one we want from time ti
				if acceptedStates[key].cum_metric < cost:
					MLkey = key
					cost = acceptedStates[key].cum_metric

		
		MLPath.append(acceptedStates[MLkey])
		# We actually only need to know where to start, so break!
		break


	# Print out r_hat
	print('After applying decoder, r_hat is then:')
	r_hat = array([])
	done = False
	qi = MLPath[0]
	for ti in range(t,-1,-1):
		key = getKey(ti,qi.prev,qi.this)

		if (ti - 1) >= 0:
			prevKey = getKey(ti-1,qs[qi.prev].prev[0],qi.prev)
			ii = 0
			if prevKey not in acceptedStates:
				prevKey = getKey(ti-1,qs[qi.prev].prev[1],qi.prev)
				ii = 1

			# Update qi for t-1
			qi = acceptedStates[prevKey]

			if qs[acceptedStates[prevKey].this].next[0] == qs[acceptedStates[key].this].this:
				bits = array([ int(x) for x in list(qs[qi.this].output[0])])
			else:
				bits = array([ int(x) for x in list(qs[qi.this].output[1])])

			# Add the bits we just found to r_hat
			r_hat = concatenate((r_hat,flip(bits,axis=0)),axis=0)

	# Flip r_hat to be the correct direction
	r_hat = flip(r_hat,axis=0)
	print(r_hat)
	print('And we found bit errors at these locations:')
	print(absolute(r - r_hat))
