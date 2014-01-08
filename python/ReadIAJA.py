__author__ = 'jdhughes'

import numpy as np
from matplotlib.pyplot import figure, show

class CrsData:
    def __init__(self, **kwargs):
        if kwargs.has_key('file'):
            self.file = kwargs['file']
        else:
            self.file = None
        self.neq = 0
        self.nnz = 0
        self.ia = []
        self.ja = []
        self.a = []
        self.rhs = []
        if self.file is not None:
            self.get_iaja(**kwargs)

    def get_iaja(self, **kwargs):
        if kwargs.has_key('file'):
            self.file = kwargs['file']
            self.neq = 0
            self.nnz = 0
            self.ia = []
            self.ja = []
        if self.file is None:
            return
        get_a = False
        if kwargs.has_key('get_a'):
            get_a = kwargs['get_a']
        get_rhs = False
        if kwargs.has_key('get_rhs'):
            get_rhs = kwargs['get_rhs']
        #--open the file
        f = open(self.file, 'r')
        #--read the data
        comment = '#'
        #--read size
        line = '#  '
        while line[0] == comment:
            line = f.readline()
            if line[0] is '#':
                continue
            t = line.split()
            self.neq, self.nnz = int(t[0]), int(t[1])
        #--dimension ia and ja
        self.ia = np.zeros(self.neq+1, np.int)
        self.ja = np.zeros(self.nnz, np.int)
        #--read ia
        line = '#  '
        while line[0] == comment:
            line = f.readline()
            if line[0] is comment:
                continue
            t = line.split()
            for idx, tv in enumerate(t):
                self.ia[idx] = int(tv)
            self.ia[-1] = self.nnz + 1
        #--read ja
        line = '#  '
        while line[0] == comment:
            line = f.readline()
            if line[0] is comment:
                continue
            t = line.split()
            for idx, tv in enumerate(t):
                self.ja[idx] = int(tv)
        #--read a vector
        if get_a is True:
            self.a = get_vector(f, self.nnz)
        #--read rhs vector
        if get_rhs is True:
            self.rhs = get_vector(f, self.neq)
        #--close the file
        f.close()

    def get_crssize(self):
        return self.neq, self.nnz

    def get_ia(self):
        return self.ia

    def get_ja(self):
        return self.ja

    def get_a(self):
        return self.a

    def get_rhs(self):
        return self.rhs

class ModelConnections:
    def __init__(self, **kwargs):
        if kwargs.has_key('file'):
            self.file = kwargs['file']
        else:
            self.file = None
        self.nconn = 0
        self.m = []
        self.n = []
        self.conductance = []
        if self.file is not None:
            self.get_model_connections(**kwargs)

    def get_model_connections(self, **kwargs):
        if kwargs.has_key('file'):
            self.file = kwargs['file']
            self.nconn = 0
            self.m = []
            self.n = []
        if self.file is None:
            return
        get_conductance = False
        if kwargs.has_key('get_conductance'):
            get_conductance = kwargs['get_conductance']
        #--open the file
        f = open(self.file, 'r')
        #--read the data
        comment = '#'
        #--read size
        line = '#  '
        while line[0] == comment:
            line = f.readline()
            if line[0] is comment:
                continue
            t = line.split()
            self.nconn = int(t[0])
        #--dimension m and n
        self.m = np.zeros(self.nconn, np.int)
        self.n = np.zeros(self.nconn, np.int)
        #--read first connection list
        self.get_connection_list(f)
        #--read second connection list
        self.get_connection_list(f)
        #--read conductance vector
        if get_conductance is True:
            self.conductance = get_vector(f, self.nconn)
        #--close the file
        f.close()

    def get_connection_list(self, f):
        comment = '#'
        line = '#  '
        while line[0] == comment:
            line = f.readline()
            if line[0] is comment:
                continue
            t = line.split()
            for idx in xrange(0,self.nconn):
                if int(t[0]) == 1:
                    self.m[idx] = int(t[idx+1])
                else:
                    self.n[idx] = int(t[idx+1])

    def get_connections(self):
        return self.nconn, self.m, self.n

    def get_conductance(self):
        return self.conductance


def plt_crs(crs_list):
    fig = figure()
    ax1 = fig.add_subplot(1,1,1)
    max_col = 0
    for [color, marker, markersize, offset, ia, ja] in crs_list:
        neq = ia.shape[0] - 1
        if (neq + offset) > max_col:
            max_col = neq + offset
        #--create empty dense matrix
        dense_matrix = np.zeros((neq+offset,neq+offset), np.int)
        #--determine if zero based or one based
        mincol = ja.min()
        idx_minus = 0
        if mincol > 0:
            idx_minus = 1
        #--fill dense matrix
        for i in xrange(0,neq):
            i0 = ia[i] - idx_minus
            i1 = ia[i+1] - idx_minus
            for j in xrange(i0,i1):
                jcol = ja[j]-idx_minus
                dense_matrix[i+offset, jcol+offset] = 1
        ax1.spy(dense_matrix, markersize=markersize, marker=marker, color=color)
    lines = np.arange(-0.5, max_col, 1)
    ax1.hlines(lines, -0.5, max_col-0.5)
    ax1.vlines(lines, -0.5, max_col-0.5)
    ax1.set_xlim(-0.5, max_col-0.5)
    ax1.set_ylim(max_col-0.5, -0.5)
    #--show the plot
    show()
    return

def get_vector(f, vlen):
    v = np.zeros(vlen, dtype=np.float)
    comment = '#'
    line = '#  '
    while line[0] == comment:
        line = f.readline()
        if line[0] is comment:
            continue
        t = line.split()
        for idx in xrange(0, vlen):
            v[idx] = float(t[idx])
    return v



