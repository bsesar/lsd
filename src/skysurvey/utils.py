import subprocess, os, errno
import numpy as np

class intervalset:
	""" A simple class representing a set of closed intervals.
	    Has two methods: add(i) to add a new interval i, and
	    isInside(x) to test whether each point in a NumPy
	    array x is inside the interval set.
	"""
	def __init__(self, *ivals):
		self.data = []
		for i in ivals:
			self.add(i)

	def add(self, i):
		# Add an interval to this interval set.
		assert len(i) == 2

		# Find intervals overlaping this one
		il = bisect_left(self.data, i[0])
		ir =  bisect_right(self.data, i[1])

		# Cut out all the points in the middle
		del self.data[il:ir]

		if il % 2 == 0:	# Left edge is outside an existing interval. Insert it
			self.data.insert(il, i[0])
			il += 1
		if ir %2 == 0: # Right edge is outside an existing interval. Insert it
			self.data.insert(il, i[1])

		assert len(self.data) % 2 == 0

	def isInside(self, x):
		""" Element-wise test if numpy array x is in the
		    interval set. Returns a bool NumPy array.
		    
		    Scales like O(len(Nivals)*len(Npts)) but since numpy
		    arrays are used, and for small interval sets, it's
		    usually fast enough (note that a Python-only "optimal"
		    O(log(Nivals)*len(Npts)) implementation would be many
		    times slower because of interpreter overhead)
		"""
		x = np.array(x, copy=False)
		ivals = np.array(self.data, copy=False)

		in_ = np.zeros(len(x), dtype=bool)
		for i in xrange(0, len(ivals), 2):
			in_ |= (ivals[i] <= x) & (x <= ivals[i+1])

		return in_

def xy_to_cell_id(x, y, t = None):
	if t is None:
		t = 0., 0., 1.
	t, t0, dt = t

	ct = astype((t - t0) / dt, np.uint64)
	ct = np.where(ct < np.uint64(0), np.uint64(0), ct)

	# construct the 32bit ID prefix from the above
	# Prefix format: 10bit x + 10bit y + 12bit time
	ix   = astype((1 + x) / 2. * 2**10, np.uint64)
	iy   = astype((1 + y) / 2. * 2**10, np.uint64)
	id   = ix << 22
	id  |= iy << 12
	id  |= ct & 0xFFF
	id <<= 32

	return id

def isiterable(x):
	try:
		iter(x)
		return True
	except TypeError:
		return False

def unpack_callable(func):
	""" Unpack a (function, function_args) tuple
	"""
	func, func_args = (func, ()) if callable(func) or func is None else (func[0], func[1:])
	return func, func_args

def gnomonic(lon, lat, clon, clat):
	from numpy import sin, cos

	phi  = np.radians(lat)
	l    = np.radians(lon)
	phi1 = np.radians(clat)
	l0   = np.radians(clon)

	cosc = sin(phi1)*sin(phi) + cos(phi1)*cos(phi)*cos(l-l0)
	x = cos(phi)*sin(l-l0) / cosc
	y = (cos(phi1)*sin(phi) - sin(phi1)*cos(phi)*cos(l-l0)) / cosc

	return (np.degrees(x), np.degrees(y))

def gc_dist(lon1, lat1, lon2, lat2):
	from numpy import sin, cos, arcsin, sqrt

	lon1 = np.radians(lon1); lat1 = np.radians(lat1)
	lon2 = np.radians(lon2); lat2 = np.radians(lat2)

	return np.degrees(2*arcsin(sqrt( (sin((lat1-lat2)*0.5))**2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)*0.5))**2 )));

_fmt_map = {
	'int8':         '%4d',
	'int16':	'%6d',
	'int32':	'%11d',
	'int64':	'%21d',
	'float32':	'%7.3f',
	'float64':	'%12.8f',
	'uint8':        '%3s',
	'uint16':	'%5s',
	'uint32':	'%10s',
	'uint64':	'%20s',
	'bool':		'%1d'
}

def get_fmt(dtype):
	#
	# Note: there's a bug with formatting long integers (they're formatted as signed), that will be fixed in numpy 1.5.1
	#       Once it's is fixed, change the format chars for uints back to 'd'
	#
	# http://projects.scipy.org/numpy/ticket/1287
	#
	if dtype.kind == 'S':
		return '%' + str(dtype.itemsize) + 's'
	if dtype == np.dtype(np.object_):
		return '%s'

	# Note: returning %s by default for unknown types
	stype = str(dtype)
	return _fmt_map[stype] if stype in _fmt_map else '%s'

def make_printf_string(row):
	fmt = ' '.join( [ get_fmt(row.dtype.fields[name][0]) for name in row.dtype.names ] )
	return fmt

def is_scalar_of_type(v, t):
	s = np.array([], v).dtype.type
	return s == t

def full_dtype(arr):
	""" Return the dtype string of the ndarray that includes
	    the array shape. Useful when merging multiple ndarray
	    columns into a single structured array table

	    See http://docs.scipy.org/doc/numpy/reference/arrays.dtypes.html for details
	"""
	shape = arr.shape[1:]
	dtype = str(shape) + str(arr.dtype) if len(shape) else str(arr.dtype)
	return dtype

def as_tuple(row):
	return tuple((row[col] for col in row.dtype.names))

def as_columns(rows, start=None, stop=None, stride=None):
	# Emulate slice syntax: only one index present
	if stop == None and stride == None:
		stop = start
		start = None
	return tuple((rows[col] for col in rows.dtype.names[slice(start, stop, stride)]))

def shell(cmd):
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	(out, err) = p.communicate();
	if p.returncode != 0:
		err = subprocess.CalledProcessError(p.returncode, cmd)
		raise err
	return out;

def mkdir_p(path):
	''' Recursively create a directory, but don't fail if it already exists. '''
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST:
			pass
		else:
			raise

def chunks(l, n):
	""" Yield successive n-sized chunks from l.
	    From http://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks-in-python
	"""
	for i in xrange(0, len(l), n):
		yield l[i:i+n]

def astype(v, t):
	""" Typecasting that works for arrays as well as scalars.
	    Note: arrays not being 1st order types in python is truly
	          annoying for scientific applications....
	"""
	if type(v) == np.ndarray:
		return v.astype(t)
	return t(v)

# extract/compact functions by David Zaslavsky from 
# http://stackoverflow.com/questions/783781/python-equivalent-of-phps-compact-and-extract
#
# -- mjuric: modification to extract to ensure variable names are legal
import inspect

legal_variable_characters = ''
for i in xrange(256):
	c = chr(i)
	legal_variable_characters = legal_variable_characters + (c if c.isalnum() else '_')

def compact(*names):
	caller = inspect.stack()[1][0] # caller of compact()
	vars = {}
	for n in names:
		if n in caller.f_locals:
			vars[n] = caller.f_locals[n]
		elif n in caller.f_globals:
			vars[n] = caller.f_globals[n]
	return vars

def extract(vars, level=1):
	caller = inspect.stack()[level][0] # caller of extract()
	for n, v in vars.items():
		n = n.translate(legal_variable_characters)
		caller.f_locals[n] = v   # NEVER DO THIS ;-)

def extract_row(row, level=1):
	caller = inspect.stack()[level][0] # caller of extract()
	for n in row.dtype.names:
		v = row[n]
		n = n.translate(legal_variable_characters)
		caller.f_locals[n] = v   # NEVER DO THIS ;-)

def xhistogram(data, bin_edges):
	""" Bin the points in 'data' into bins whose edges are given by bin_edges.
	    The output array at location i will contain the number of points pts
	    satisfying bin_edges[i-1] < pts < bin_edges[i]
	    
	    Points less than bin_edges[0] and greater than bin_edges[-1] will be
	    at indices 0 and len(bin_edges) in the output array, respectively.
	"""
	bins = np.empty(len(bin_edges)+2, dtype='f8')
	bins[0]    = -np.inf
	bins[1:-1] = bin_edges
	bins[-1]   =  np.inf
	hist, _ = np.histogram(data, bins)
	return hist
