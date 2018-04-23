## Viterbi Decoder
# Nicholas McKibben
# ECEn 770
# 2018-04-14

from numpy import *
from scipy import special as sp
import matplotlib.pyplot as plt

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


def qfunc(xlist):
	# val = 0.5 - 0.5*sp.erf(xlist/sqrt(2))
	# return(val)
	if type(xlist) is float64:
		xlist = [ xlist ]

	p = zeros(len(xlist))
	for ii in range(len(xlist)):
		if xlist[ii] < 0:
			p[ii] = 1 - 0.5*sp.erfc(-xlist[ii]/sqrt(2));
		else:
			p[ii] = 0.5*sp.erfc(xlist[ii]/sqrt(2));

	return(p)

def crossprob(N0,Ec):
	#print('N0: %f, Ec: %f' % (N0,Ec))
	val = sqrt(2*Ec/N0)
	val = qfunc(val)
	return(val)

## Convolutional Encoder
def convencode(m,g):
    # Convolve for each g[j]
    c = zeros([ g.shape[0],(m.shape[0] + g.shape[1] -1) ])
    for ii in range(0,g.shape[0]):
        c[ii,:] = mod(convolve(m,g[ii,:]),2)
    
    # Interleave the rows
    c = reshape(c,[ 1,-1 ],order='F')
    return(c[0])

def decode(r,qs,hard=True):
	# Set up occupied states
	currStates = Occupied()
	currStates.add([ 0 ]) # start in the 0th state

	# Set up accepted states array
	acceptedStates = dict()
	key = getKey(0,0,0)
	acceptedStates[key] = Accepted(0,0,0,0)

	# Set up the branches
	branches = []

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

	# De-interleave r_hat
	r_hat0 = r_hat[::2]
	r_hat1 = r_hat[1::2]

	# Now find m_hat from r_hat
	m_hat = mod(r_hat0[1::] + r_hat1[1::],2)
	
	return(r_hat,m_hat)

if __name__ == '__main__':

	## (1) Encoder
	# Write a computer program that performs the encoding operation
	# for a convolutional code with transfer matrix:
	g1 = array([ 1,0,1 ])
	g2 = array([ 1,1,1 ])

	# Test the encoder
	m = array([ 1,1,0,0,1,0,1 ])
	# m = random.randint(2,size=100000)
	c = convencode(m,array([ g1,g2 ]))
	c = c[:-2]
	print(c)

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
	# n = [ int(x < 0.5) for x in random.random(size=len(c)) ]
	# r = mod(c + n,2)
	print('Start with r:')
	print(r.astype(float))

	# Get some states
	q1 = State(0,[ 0,1 ],[ 0,2 ],[ '00','11' ])
	q2 = State(1,[ 2,3 ],[ 0,2 ],[ '01','10' ])
	q3 = State(2,[ 0,1 ],[ 1,3 ],[ '11','00' ])
	q4 = State(3,[ 2,3 ],[ 1,3 ],[ '10','01' ])
	qs = [ q1,q2,q3,q4 ]

	# Perform the decoding
	r_hat,m_hat = decode(r,qs,hard=True)

	print('After applying decoder, r_hat is then:')
	print(r_hat)
	print('And we found bit errors at these locations:')
	print(absolute(r - r_hat))

	print('Recall that m was:')
	print(m.astype(float))
	print('And we get m_hat is:')
	print(m_hat)
	print('Bit Error is %d' % sum(mod(m + m_hat,2)))


	# input()

	## Now Simulate
	k = 10000
	trace_back = 20
	rate = 1/2
	Ec = 1
	Eb = Ec/rate
	gammas = linspace(1,6,5)
	N = 100
	sigma2 = zeros(len(gammas))
	N0 = zeros(len(gammas))
	Pe = zeros(len(gammas))

	for ii in range(len(gammas)):
		# Compute N0, sigma^2
		N0[ii] = Ec/(rate*gammas[ii])
		sigma2[ii] = N0[ii]/2
		p = crossprob(N0[ii],Ec)

		nn = 0
		nbits = 0
		while nn < N:
			# Generate message and codeword. Assume random.
			m = random.randint(2,size=k)
			c = convencode(m,array([ g1,g2 ]))
			c = c[:-2]
			#print('c')
			#print(c)

			# Generate noise and recieved signal r
			n = [ int(x < p) for x in random.random(size=len(c)) ]
			#print('n, p = %f' % p)
			#print(n)
			r = mod(c + n,2)
			#print('r')

			# Add some bits
			nbits += k

			# Decode the recieved codeword and get m_hat
			r_hat,m_hat = decode(r,qs,hard=True)
			#print(mod(m_hat + m,2))

			# Accumulate error
			nn += sum(mod(m[0:trace_back] + m_hat[0:trace_back],2))
			#print(mod(m+m_hat,2))
			# print(nn)
			# input()

		Pe[ii] = nn/nbits
		print('gamma = %f is done with Pe = %f!' % (gammas[ii],Pe[ii]))

	# Get Theoretical uncoded BPSK
	BPSK = qfunc(sqrt(2*Eb/N0[:]))

	# Plot it
	plt.semilogy(10*log10(Eb/N0),Pe)
	plt.semilogy(10*log10(Eb/N0),BPSK)
	plt.grid(True)
	plt.title('Bit Error')
	plt.xlabel('E_b/N_0')
	plt.ylabel('P_e')
	plt.show()