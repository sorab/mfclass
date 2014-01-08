#todo list:
#   1.  DONE--Merge submodel.py and mf2015.py into one main program
#       that works for either implicit or explicit cross
#   2.  DONE--get submodel.py reading from external files
#   3.  DONE--Get picard concepts fleshed out.  Remove models from
#       main and make solution.solve
#   4.  Put in a time step loop
#   5.  Convert to FORTRAN


import numpy as np
from ReadIAJA import *
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve, cg

class Sparse(object):
    def __init__(self, nrow, rowmaxnnz):
        self.nrow = nrow
        self.row = np.empty( (nrow), dtype=object)
        self.nnz = 0
        for n in xrange(self.nrow):
            self.row[n] = []
        return

    def addconnection(self, n, m):
        self.row[n].append(m)
        self.nnz += 1
        return

    def getiaja(self, ia, ja):
        ipos = 0
        for n in xrange(self.nrow):
            ia[n] = ipos
            columnlist = self.row[n]
            if n in columnlist: #diagonal
                ja[ipos] = n
                columnlist.remove(n)
                ipos += 1
            columnlist.sort()
            for j in columnlist:
                ja[ipos] = j
                ipos += 1
        ia[self.nrow] = ipos
        return ia, ja


class Solution(object):
    def __init__(self, mxiter, xtol):
        #create an empty solution
        self.models = []
        self.crosses = []
        self.mxiter = mxiter
        self.xtol = xtol
        self.neq = 0
        self.nja = 0
        return

    def add_model(self, model):
        """
        Add a model to this solution

        """
        self.models.append(model)
        return

    def add_cross(self, cross):
        """
        Add a cross to this solution

        """
        self.crosses.append(cross)
        return

    def initialize(self):

        #go through each model and assign offsets and solution
        for m in self.models:
            m.offset = self.neq
            self.neq += m.neq
            m.solution = self

        #allocate arrays
        self.ia = np.empty( (self.neq + 1), dtype=np.int)
        self.x = np.zeros( (self.neq), dtype=np.float)
        self.rhs = np.empty( (self.neq), dtype=np.float)

        #placeholders for arrays to be filled
        self.nja = None
        self.ja = None
        self.amat = None

        #go through each model and point m.x and m.rhs
        for m in self.models:
            m.x = self.x[m.offset:m.offset+m.neq]
            m.rhs = self.rhs[m.offset:m.offset+m.neq]

        #initialize and store the sparse object
        self.sparse = Sparse(self.neq, 0)

        return

    @staticmethod
    def create(filename):
        f = open(filename, 'r')
        lines = f.readlines()
        for line in lines:
            if line.strip() == '':
                continue
            varname, varval = line.strip().split()
            if 'mxiter' in varname:
                mxiter = int(varval)
            elif 'xtol' in varname:
                xtol = float(varval)
        f.close()
        s = Solution(mxiter, xtol)
        return s

    def reset(self):
        """
        Set amat and rhs to zero.

        """
        self.amat[:] = 0.
        self.rhs[:] = 0.
        return

    def solve(self):
        icnvg = 0
        for kiter in xrange(self.mxiter):
            print 'SOLVING ITERATION: ', kiter + 1
            self.reset()
            for m in self.models:
                m.fill()
            for m in self.models:
                m.fmcalc()
            for c in self.crosses:
                c.fill()
            icnvg = self.linsolve()
            if icnvg == 1:
                break
        return

    def linsolve(self):
        icnvg = 0
        xsaved = self.x.copy()
        acsr = csr_matrix( (self.amat, self.ja, self.ia), shape=(self.neq,
                                                                 self.neq) )
        info = 0
        #can't use spsolve here because it sorts ja (diagonal not in front)
        #self.x[:] = spsolve(acsr, self.rhs)
        self.x[:], info = cg(acsr, self.rhs, maxiter=2000)
        if abs(self.x - xsaved).max() < self.xtol:
            icnvg = 1
        if info < 0:
            print 'illegal input or breakdown in linear solver...'
            raise Exception()
        return icnvg

    def connect(self):
        
        #add the internal model connections to self.sparse
        #MOVE TO A NEW MODEL.CONNECT METHOD
        for m in self.models:
            for i in xrange(m.neq):
                for jj in xrange(m.ia[i], m.ia[i+1]):
                    j = m.ja[jj]
                    self.sparse.addconnection(i + m.offset, j + m.offset)
            
        #add the cross terms to self.sparse
        #MOVE TO A NEW CROSS.CONNECT METHOD
        for c in self.crosses:
            if not c.implicit:
                continue
            m1 = c.m1
            m2 = c.m2
            ncross = c.ncross
            nodem1 = c.nodem1
            nodem2 = c.nodem2
            for n in xrange(ncross):
                i = c.nodem1[n]
                j = c.nodem2[n]
                self.sparse.addconnection(i + m1.offset, j + m2.offset)
                self.sparse.addconnection(j + m2.offset, i + m1.offset)

        #create the global ja array
        #THIS STAYS HERE AS IT IS WORKING ON ITSELF
        self.nja = self.sparse.nnz
        self.ja = np.empty( (self.nja), dtype=np.int)
        self.amat = np.empty( (self.nja), dtype=np.float)
        self.ia, self.ja = self.sparse.getiaja(self.ia, self.ja)
        
        #create the arrays for mapping model ja positions to global solution
        #MOVE TO A NEW MODEL.IDXGLOMAP METHOD
        for m in self.models:
            ipos = 0
            for n in xrange(m.neq):
                iabeg = self.ia[n + m.offset]
                iaend = self.ia[n + m.offset + 1]
                for jj in xrange(iabeg, iaend):
                    j = abs(self.ja[jj])
                    if m.offset <= j < m.offset + m.neq:
                        m.idxglo[ipos] = jj
                        ipos += 1
                        
        #create the arrays for mapping cross connections to global solution
        #MOVE TO A NEW CROSS.IDXGLOMAP METHOD
        for c in self.crosses:
            if not c.implicit:
                continue
            m1 = c.m1
            m2 = c.m2
            idxglo = c.idxglo
            idxsymglo = c.idxsymglo
            ncross = c.ncross
            nodem1 = c.nodem1
            nodem2 = c.nodem2
            for n in xrange(ncross):
                iglo = nodem1[n] + m1.offset
                jglo = nodem2[n] + m2.offset
                #find the jglobal value in row iglo and store in idxglo
                for ipos in xrange(self.ia[iglo], self.ia[iglo + 1]):
                    if jglo == self.ja[ipos]:
                        idxglo[n] = ipos
                        break
                #find the symmetric counterpart and store in idxsymglo
                for ipos in xrange(self.ia[jglo], self.ia[jglo + 1]):
                    if iglo == self.ja[ipos]:
                        idxsymglo[n] = ipos
                        break
               
        return


class Model(object):
    def __init__(self, name, neq, ia, ja, cond1, rhs1):
        self.name = name
        self.neq = neq
        self.nja = ja.shape[0]
        self.ia = ia
        self.ja = ja
        self.cond = cond1
        self.rhs = rhs1
        self.idxglo = -1 * np.ones( (self.nja), dtype=np.int)
        self.solution = None
        self.offset = 0
        self.packages = []
        self.models = []
        self.crosses = []
        self.solution = None
        return

    @staticmethod
    def create(filename):
        f = open(filename, 'r')
        modelname = filename.strip().split('.')[0]
        for line in f:
            linesplit = line.strip().split()
            ftype, fname = linesplit[0:2]
            if 'dis' in ftype:
                m1_iaja = CrsData(file=fname, get_a=True,get_rhs=True)
                neq1, nnz1 = m1_iaja.get_crssize()
                ia1 = m1_iaja.get_ia()
                ja1 = m1_iaja.get_ja()
                cond1 = m1_iaja.get_a()
                rhs1 = m1_iaja.get_rhs()
                ia1 = ia1 - 1
                ja1 = ja1 - 1
                m = Model(modelname, neq1, ia1, ja1, cond1, rhs1)
            elif 'wel' in ftype:
                p = WEL.create(fname, m)
            elif 'ghb' in ftype:
                p = GHB.create(fname, m)
        f.close()
        return m

    def fill(self):
        print 'FMFILL...', self.name
        amatglo = self.solution.amat
        rhsglo = self.solution.rhs
        #fill amat with conductance from this model
        for n in xrange(self.nja):
            amatglo[self.idxglo[n]] = self.cond[n]
        #fill global rhs with conductance from this model
        #not needed because self.rhs points to rhsglo
        #for n in xrange(self.neq):
        #    rhsglo[n + self.offset] += self.rhs[n]
        return

    def fmcalc(self):
        print 'FMCALC...', self.name
        for package in self.packages:
            package.fmcalc()
            package.fmfill()
        return


class Package(object):
    def __init__(self, name, model, nbound, nodelist):
        self.name = name
        self.model = model
        model.packages.append(self)
        self.nbound = nbound
        self.nodelist = nodelist
        self.hcof = np.empty((nodelist.shape[0]), dtype=np.float)
        self.rhs = np.empty((nodelist.shape[0]), dtype=np.float)
        return
        
    def fmfill(self):
        print 'FM PACKAGE FILL...', self.name
        for i in xrange(self.nbound):
            node = self.nodelist[i]
            self.model.rhs[node] += self.rhs[i]
            nodeglo = node + self.model.offset
            #add hcoff to diagon of parent models solution
            ipos = self.model.solution.ia[nodeglo]
            self.model.solution.amat[ipos] += self.hcof[i]
        return


class GHB(Package):
    def __init__(self, name, model, nbound, nodelist, stage, cond):
        super(GHB, self).__init__(name, model, nbound, nodelist)
        self.cond = cond
        self.stage = stage
        return

    @staticmethod
    def create(filename, model):
        packagename = filename  #.strip().split('.')[0]
        f = open(filename, 'r')
        line = f.readline()
        mxbnd = int(line)
        line = f.readline()
        nbound = int(line)
        nodelist = np.empty((nbound),dtype=np.int)
        stage = np.empty((nbound),dtype=np.float)
        cond = np.empty((nbound),dtype=np.float)
        for i in xrange(nbound):
            line = f.readline()
            linelist = line.strip().split()
            nodelist[i] = np.int(linelist[0]) - 1
            stage[i] = np.float(linelist[1])
            cond[i] = np.float(linelist[2])
        f.close()
        p = GHB(packagename, model, nbound, nodelist, stage, cond)
        return p

    def fmcalc(self):
        print 'FM PACKAGE CALC...', self.name
        for i in xrange(self.nbound):
            self.hcof[i] = -self.cond[i]
            self.rhs[i] = -self.cond[i] * self.stage[i]
        return


class WEL(Package):
    def __init__(self, name, model, nbound, nodelist, q):
        super(WEL, self).__init__(name, model, nbound, nodelist)
        self.q = q
        return

    @staticmethod
    def create(filename, model):
        packagename = filename  #.strip().split('.')[0]
        f = open(filename, 'r')
        line = f.readline()
        mxbnd = int(line)
        line = f.readline()
        nbound = int(line)
        nodelist = np.empty((nbound),dtype=np.int)
        q = np.empty((nbound),dtype=np.float)
        for i in xrange(nbound):
            line = f.readline()
            linelist = line.strip().split()
            nodelist[i] = np.int(linelist[0]) - 1
            q[i] = np.float(linelist[1])
        f.close()
        p = WEL(packagename, model, nbound, nodelist, q)
        return p

    def fmcalc(self):
        print 'FM PACKAGE CALC...', self.name
        for i in xrange(self.nbound):
            self.hcof[i] = 0
            self.rhs[i] = -self.q[i]
        return


class Cross(object):
    def __init__(self, name, m1, m2, ncross, nodem1, nodem2, cond):
        self.name = name
        self.m1 = m1
        self.m2 = m2
        self.implicit = True
        self.ncross = ncross
        self.nodem1 = nodem1
        self.nodem2 = nodem2
        self.cond = cond
        self.idxglo = -1 * np.ones( (self.ncross), dtype=np.int)
        self.idxsymglo = -1 * np.ones( (self.ncross), dtype=np.int)
        self.solution = None
        return

    def initialize(self):
        """
        Determine if cross is implicit or explicit based on whether or not
        m1 and m2 are part of the same solution
        """
        if self.m1.solution != self.m2.solution:
            self.implicit = False
            self.m1.solution.add_cross(self)
            self.m2.solution.add_cross(self)

            cond = np.zeros((self.ncross), dtype=np.float)

            name = self.m2.name + ' to ' + self.m1.name
            stage1 = np.zeros((self.ncross), dtype=np.float)
            self.m1xrsp = GHB(name, self.m1, self.ncross, self.nodem1, stage1,
                              self.cond)

            name = self.m1.name + ' to ' + self.m2.name
            stage2 = np.zeros((self.ncross), dtype=np.float)
            self.m2xrsp = GHB(name, self.m2, self.ncross, self.nodem2, stage2,
                              cond)
        else:
            #m1 and m2 part of same solution so add self to the solution
            #and point self.solution to m1.solution
            self.m1.solution.add_cross(self)
            self.solution = self.m1.solution
        return

    @staticmethod
    def create(filename, m1, m2):
        name = filename.strip().split('.')[0]
        m1_m2_conn = ModelConnections(file=filename, get_conductance=True)
        ncross, xrsnodem1, xrsnodem2 = m1_m2_conn.get_connections()
        xrscond = m1_m2_conn.get_conductance()
        xrsnodem1 -= 1
        xrsnodem2 -= 1
        c = Cross(name, m1, m2, ncross, xrsnodem1, xrsnodem2, xrscond)
        return c

    def fill(self):
        if self.implicit:
            self.fill_implicit(self.cond)
        else:
            self.fill_explicit(self.cond)
        return

    def fill_explicit(self, cond):
        #substitute x into stage for model 1
        for i in xrange(self.ncross):
            n1 = self.nodem1[i]
            n2 = self.nodem2[i]
            self.m1xrsp.stage[i] = self.m2.x[n2]
            self.m2xrsp.stage[i] = self.m1.x[n1]
            self.m1xrsp.cond[i] = cond[i]
        return
        
    def fill_implicit(self, cond):
        amatglo = self.solution.amat
        for n in xrange(self.ncross):
            amatglo[self.idxglo[n]] = cond[n]
            amatglo[self.idxsymglo[n]] = cond[n]  #symmetric counterpart
            nodem1glo = self.nodem1[n] + self.m1.offset
            nodem2glo = self.nodem2[n] + self.m2.offset
            idiaglo = self.solution.ia[nodem1glo]
            amatglo[idiaglo] -= cond[n]
            idiaglo = self.solution.ia[nodem2glo]
            amatglo[idiaglo] -= cond[n]
        return

def namefile(fname):
    """
    Read the name file and create instances
    """
    models = []
    crosses = []
    solutions = []
    solution_groups = []

    f = open('simulation.nam', 'r')
    while True:
        try:
            line = f.next()
        except:
            break
        if '#' in line[0]:
            continue
        elif line.strip() == '':
            continue
        elif 'model' in line[0:5]:
            ftype, filename, modelid = line.strip().split()
            m = Model.create(filename)
            models.append(m)
        elif 'xrs' in line[0:3]:
            ftype, filename, m1id, m2id = line.strip().split()
            m1id = int(m1id)
            m2id = int(m2id)
            c = Cross.create(filename, models[m1id - 1], models[m2id - 1])
            crosses.append(c)
        elif 'sms' in line[0:3]:
            ftype, filename = line.strip().split()[0:2]
            s = Solution.create(filename)
            solutions.append(s)
        elif 'tdis' in line[0:4]:
            print 'tdis not implemented yet'
        elif 'begin solution_group' in line[0:20]:
            solutiongrouplist = []
            line = f.next()
            varname, varval = line.strip().split()[0:2]
            if 'mxiter' in varname:
                mxiter = int(varval)
            line = f.next()
            varname, varval = line.strip().split()[0:2]
            if 'nsolutions' in varname:
                nsolutions = int(varval)
            for n in xrange(nsolutions):
                line = f.next()
                varname, varval = line.strip().split()[0:2]
                if 'sms' in varname:
                    isolnid = int(varval) - 1
                    s = solutions[isolnid]
                line = f.next()
                varname, varval = line.strip().split()[0:2]
                if 'nmodels' in varname:
                    nmodels = int(varval)
                for imodels in xrange(nmodels):
                    modelgrouplist = []
                    line = f.next()
                    varname, varval = line.strip().split()[0:2]
                    if 'model' in varname:
                        mid = int(varval) - 1
                        m = models[mid]
                        modelgrouplist.append(m)
                        s.add_model(m)
                solutiongrouplist.append([s, modelgrouplist])
            sg = SolutionGroup(mxiter, nsolutions, solutiongrouplist)
            solution_groups.append(sg)
    f.close()

    print 'nsolution_groups', len(solution_groups)
    print 'nsolutions: ', len(solutions)
    print 'nmodels: ', len(models)
    print 'ncrosses: ', len(crosses)

    return models, crosses, solutions, solution_groups

class SolutionGroup(object):
    def __init__(self, mxiter, nsolutions, solutions):
        """
        A group of solutions that are coupled
        solutions = [ [s1, [m1, m2]], [s2, [m3, m4]], ...]

        """
        self.mxiter = mxiter
        self.nsolutions = nsolutions
        self.solutions = solutions
        return

    def solve(self):
        """
        Solve the solution group

        """
        for kiter in xrange(self.mxiter):
            for s, modelgroup in self.solutions:
                print 'SOLUTION GROUP ITERATION', kiter + 1
                s.solve()
        return

def main():
    models, crosses, solutions, solution_groups = namefile('simulation.nam')

    for s in solutions:
        s.initialize()

    for c in crosses:
        c.initialize()

    for s in solutions:
        s.connect()

    for sg in solution_groups:
        sg.solve()

    s = solutions[0]
    print 'amat1: ', s.amat
    print 'ja: ', s.ja
    print 'rhs: ', s.rhs
    print 'answer: ', s.x
    #print 'model1 x: ', models[1-1].x
    #print 'model2 x: ', models[2-1].x
    return


if __name__ == "__main__":
    main()

