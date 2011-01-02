#!/usr/bin/env python

import numpy as np
import sys
import mrp2p.peer

class MyClass(object):
	nada = "aaa"
	pass

def mapper(values):
	for v in values:
		my = MyClass()
		my.ct = 1
		yield v % 512, (my, np.arange(1024))

def reducer1(kv):
	k, v = kv
	key = 'even' if k % 2 == 0 else 'odd'
	yield key, sum(my.ct for my, _ in v)

def reducer(kv):
	k, v = kv
	yield k, sum(v)

if __name__ == "__main__":
	v = np.arange(4*2**10)

	# Parallel
	pool = mrp2p.peer.Pool('peers')
	res1 = []
	for res in pool.map_reduce_chain(v, [mapper, reducer1, reducer]):
		print >>sys.stderr, "Result: ", res
		res1.append(res)
	res1 = sorted(res1)
	exit()

	# Classic
	res2 = sorted(list(mapper(v)))

	if res1 == res2:
		print "Result OK."
	else:
		print res1
		print res2
		print "Check FAILED"
