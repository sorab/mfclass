#cdl, 1/21/2014
#simple test program to read fortran source files.
#it uses the crackfortran module


import os
import numpy.f2py.crackfortran as cf
srcpath = '../src'
filelist = os.listdir(srcpath)
print filelist
filelist2=[]
for fname in filelist:
    filelist2.append(srcpath + '/' + fname)

for fname in filelist2:
    print fname

postlist = cf.crackfortran(filelist2[:-1])

def dimlst2str(dimlst):
    s = '('
    for dm in dimlst:
        s += dm + ','
    s = s[:-1] + ')'
    if 'njas' in s: print dimlst
    return s
    
for b in postlist:
    print 'Blockname: ', b['name']
    if b.has_key('vars'):
        for varname in b['vars'].keys():
            d = b['vars'][varname]
            s = ''
            s += varname + ' :: '
            if d.has_key('typespec'):
                s += d['typespec'] + ' '
            if d.has_key('dimension'):
                s += dimlst2str(d['dimension'])
            print '  Variable: ', s
