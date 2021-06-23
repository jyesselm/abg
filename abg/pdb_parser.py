"""
pdblib, version 1.2, 12-09-2014
Written by Yi Xue
Copyright: Yi Xue and Skrynnikov's group @Purdue University

A simple python package to manipulate pdb files
pdb.base has not dependence of numpy, and thus does not support advanced
operations such as concatenating two chains, rotating a group of atoms, etc

"""

#import pdb  #for debug
from sys import stdout,exit
from operator import add
from math import sqrt
from functools import reduce

from abg.base import cl,pager,divide,partition
from abg.bio import aa_abbr,resabbr


#====== Atom ===================================================================
class Atom:
    def __init__(self):
        self.atid = 0
        self.loc = ' '
        self.r = None
        self.charge = None
        self.elem = ''
        self.sgid = '    '
        self.chid = ' '
        self.oc = 1.0
        self.bf = 0.0

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def gettype(self):
        return (self.name[1] if self.name[0].isdigit() else self.name[0])

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def show(self, disp=True):
        output = []
        fmtstr = '%-4s %5d %s'
        attrs = [self.name, self.atid, self.loc]
        if self.r:
            fmtstr += ' '+'%8.3f'*3
            attrs += [self.r[0], self.r[1], self.r[2]]
        else:
            fmtstr += ' '+'   *.***'*3
        if self.charge:
            fmtstr += '  %8.3f'
            attrs += [self.charge]
        else:
            fmtstr += '     *.***'
        output.append((fmtstr+'\n')%tuple(attrs))

        if disp:
            stdout.writelines(output)
        else:
            return output

#====== Residue ================================================================
class Residue:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self):
        self.name = ''
        self.resi = 0
        self.atoms = []

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # disp=True: display on screen
    def show(self, fmt='S', disp=True):
        output = []
        try:
            uid = self.uid
        except AttributeError:
            uid = ''
        output.append('%s%03d %-4s %-4s%s\n'%(cl.y, self.resi, self.name, uid,
                      cl.n))

        if fmt.upper() == 'S':
            atstrs = [(at.name+(4-len(at.name))*' ' if at.r else
                       cl.lr+at.name+(4-len(at.name))*' '+cl.n)
                      for at in self.atoms]
            output += ' '.join(atstrs) + '\n'
        else:
            for at in self.atoms:
                output += at.show(disp=False)

        if disp:
            stdout.writelines(output)
        else:
            return output

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Given atom name, return atom object; On failure, return None
    def getat(self, name):
        ats = list(filter(lambda x: x.name==name, self.atoms))
        if ats:
            return ats[0]
        else:
            return None


#====== Segment ================================================================
class Segment:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self, reses=None):
        self.sgid = '    '
        self.chid = ' '
        self.nter = False
        self.cter = False

        if reses is None:
            self.reses = []
        else:
            self.reses = reses

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def show(self, fmt='S', disp=True):
        output = ['']
        output.append(cl.g + 'Segment ' + self.sgid +
                      ', # of residue: %d'%len(self.reses) + cl.n + '\n')
        output.append(cl.g + '-'*80 + cl.n + '\n')
        if not self.reses:
            output.append('No residue exists!\n')
        else:
            for res in self.reses:
                output += res.show(fmt, disp=False)

        if disp:
            pager(output)
        else:
            return output


    def seq(self, fmt='s', dsp=1):    
        #-----------------------------------------------------------------------
        def abbr(resname):
            try:
                name = resabbr[resname]
            except KeyError:
                name = 'X'
            return name
        #-----------------------------------------------------------------------

        resns = map(lambda x: x.name, self.reses)
        if fmt.upper() == 'S':
            return ''.join(map(abbr, resns))
        else:
            return resns

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getats(self):
        return reduce(add, map(lambda x: x.atoms, self.reses))

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getindex(self, resi):
        try:
            return map(lambda x: x.resi, self.reses).index(resi)
        except ValueError:
            return None

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def renumber(self, iat=1, ires=None):
        if iat is not None:
            ats = list(filter(lambda x: x.r, self.getats()))
            for i,at in enumerate(ats):
                at.atid = i+iat
        if ires is not None:
            reses = self.reses
            for i,res in enumerate(reses):
                res.resi = i+ires
                for at in res.atoms:
                    at.resi = res.resi


#====== Molecule (or Model) ====================================================
class Mol:
    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self, inp=None):
        self.fmt = 'plain'
        self.mdid = 0
        self.top = None
        self.segs = []
        if isinstance(inp, str):
            # inp is pdb file name
            self.read(inp)
        elif isinstance(inp, Segment):
            # inp is segment, then build a Mol based on this single seg
            self.segs = [inp]

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def read(self, inp):
        nt_atom = set(['HT3','HN3','H3','3HN','3H'])
        ct_atom = set(['OT2','OXT','OB'])
        fmt = 'pdb'
        if isinstance(inp, str):
            if inp[-4:] == '.pqr':
                fmt = 'pqr'
            lines = open(inp).readlines()
            kwds=('ATOM','HETATM','TER', 'MODEL','ENDMDL','END')
            ##'Q' is pseudo atoms in molmol format pdb files
            lines = list(filter(lambda x:x[:6].rstrip() in kwds and x[13:14]!='Q',
                           lines))
            ##lgrps: lines of groups; lmds: text lines of models
            lgrps = partition(lines, lambda x: x[:6].rstrip()=='END')
            lmds = reduce(add, [partition(lgrp, lambda x: x[:5]=='MODEL',
                                include='header') for lgrp in lgrps])
            inp = lmds[0]
            if inp[0][:5]=='MODEL':
                self.mdid = int(inp[0].split()[1])
                inp = inp[1:]
            if inp[-1][:6]=='ENDMDL':
                inp = inp[:-1]

        # reading a prot/model
        self.segs = []
        ##lgrps: lines of groups
        lgrps = partition(inp, lambda x: x[:3]=='TER')
        lsegs = reduce(add,
                       [divide(lgrp, lambda x: x[72:76]) for lgrp in lgrps])
        for lseg in lsegs:
            # Building a segment
            seg = Segment()
            lreses = divide(lseg, lambda x: x[22:27])
            for lres in lreses:
                # Building a residue
                res = Residue()
                if fmt == 'pqr':
                    res.atoms = [readatom(line, fmt='pqr') for line in lres]
                else:
                    res.atoms = [readatom(line) for line in lres]
                res.name = res.atoms[0].resn
                res.resi = res.atoms[0].resi
                res.icode = res.atoms[0].icode
                seg.reses.append(res)
            atns = set(map(lambda x: x.name,seg.reses[0].atoms))
            # assume nter if is_aa and (no H or containing nt_atom)
            if seg.reses[0].name in aa_abbr and \
               ('H' not in map(Atom.gettype, seg.reses[0].atoms) \
               or (nt_atom & atns)):
                seg.nter = True
            atns = set(map(lambda x: x.name,seg.reses[-1].atoms))
            # assume cter if is_aa and (no H or containing ct_atom)
            if seg.reses[-1].name in aa_abbr and \
               ('H' not in map(Atom.gettype, seg.reses[-1].atoms) \
               or (ct_atom & atns)):
                seg.cter = True
            seg.sgid = seg.reses[0].atoms[0].sgid
            seg.chid = seg.reses[0].atoms[0].chid
            self.segs.append(seg)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # update atom entry for these fields: resn, resi, icode, chid, sgid
    def atsync(self):
        for seg in self.segs:
            for res in seg.reses:
                for at in res.atoms:
                    at.resn = res.name
                    at.resi = res.resi
                    at.icode = res.icode
                    at.chid = seg.chid
                    at.sgid = seg.sgid

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def write(self, out, header=['REMARK CREATED BY PDBLIB\n'], sync=True):
        if isinstance(out, str):
            fout = open(out, 'wt')
        else:
            fout = out

        if header:
            fout.writelines(header)

        if sync:
            self.atsync()

        for seg in self.segs:
            for res in seg.reses:
                for at in res.atoms:
                    if at.r:
                        fout.write(writeatom(at))
            fout.write('TER\n')

        if isinstance(out, str):
            fout.write('END\n')
            fout.close()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def show(self, fmt='S', disp=True):
        output = []
        output.append(('Model %d,'%self.mdid if self.mdid else 'Mol,') +
                       ' # of segment: %d\n'%len(self.segs))
        output.append('='*80 + '\n')
        for seg in self.segs:
            output += seg.show(fmt, disp=False)

        if disp:
            pager(output)
        else:
            return output

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def renumber(self, iat=1, ires=None):
        if iat is not None:
            ats = list(filter(lambda x: x.r, self.getats()))
            for i,at in enumerate(ats):
                at.atid = i+iat
        if ires is not None:
            reses = self.getreses()
            for i,res in enumerate(reses):
                res.resi = i+ires
                for at in res.atoms:
                    at.resi = res.resi

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getats(self):
        return reduce (add, map(Segment.getats, self.segs))

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getreses(self):
        return reduce (add, map(lambda x: x.reses, self.segs))


#====== Pdb ====================================================================
class Pdb:
    "A simple PDB python lib written by Yi Xue @ PULSe, Purdue Unversity"

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def __init__(self, fn=None):
        self.mds = []
        if fn:
            self.read(fn)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def read(self, fn):
        lines = open(fn).readlines()
        kwds=('ATOM','HETATM','TER','MODEL','ENDMDL','END')
        ##'Q' is pseudo atoms in molmol format pdb files
        lines = list(filter(lambda x:x[:6].rstrip() in kwds and x[13:14]!='Q', lines))

        self.mds = []
        ##lgrps: lines of groups; lmds: text lines of models
        lgrps = partition(lines, lambda x: x[:6].rstrip()=='END')
        lmds = reduce(add, [partition(lgrp, lambda x: x[:5]=='MODEL',
                            include='header') for lgrp in lgrps])
#       pdb.set_trace() #<============ for debug
        for lmd in lmds:
            md = Mol()
            if lmd[0][:5]=='MODEL':
                md.mdid = int(lmd[0].split()[1])
                lmd = lmd[1:]
            if lmd[-1][:6]=='ENDMDL':
                lmd = lmd[:-1]
            md.read(lmd)
            self.mds.append(md)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # update atom entry for these fields: resn, resi, chid, sgid
    def atsync(self):
        for md in self.mds:
            md.atsync()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # header=1: Print a header line
    def write(self, fn, header=['REMARK CREATED BY PDBLIB\n'], sync=True):
        fout = open(fn, 'wt')

        if header:
            fout.writelines(header)

        if sync:
            self.atsync()

        for md in self.mds:
            if len(self.mds) > 1:
                fout.write('MODEL%9d\n'%md.mdid)
            md.write(fout, header=[], sync=False)
            if len(self.mds) > 1:
                fout.write('ENDMDL\n')

        fout.write('END   \n')
        fout.close()

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def show(self, fmt='S'):
        output = []
        output.append(cl.y + 'PDB file, # of model: %d'%len(self.mds) +
                      cl.n + '\n')
        output.append(cl.y + '+'*80 + cl.n + '\n')
        for md in self.mds:
            output += md.show(fmt, disp=False)
        pager(output)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def format(self, fmt='auto'):
        if self.mds:
            self.mds[0].format(fmt)
        for md in self.mds[1:]:
            md.top = self.mds[0].top
            md.format(fmt)


    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def renumber(self, iat=1, ires=None):
        for md in self.mds:
            md.renumber(iat, ires)

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getats(self):
        return reduce (add, map(Mol.getats, self.mds))

    # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    def getreses(self):
        return reduce (add, map(Mol.getreses, self.mds))


        
#====== Public Functions =======================================================
# ==============================================================================
def readatom(line, fmt='pdb'):
    """
    Pdb format:
    1-6:   "ATOM  ";     7-11: ATOM ID
    13-16: Atom Name;    17: Location indicator
    18-20: resname;      22: chain identifier
    23-26: resnum;       27: icode
    31-38: X,real(8.3)   39-46: Y,real(8.3)
    47-54: Z,real(8.3)   55-60: Occupancy
    61-66: TempFactor    73-76: segID
    77-78: element       79-80: Charge
    """
    if fmt == 'pqr':
        try:
            atom = Atom()
            fds = line.split()
            atom.atid = int(fds[1])
            atom.name = fds[2]
            atom.resn = fds[3]
            atom.resi = int(fds[4])
            atom.r = tuple(map(float, fds[5:8]))
            atom.charge = float(fds[8])
            atom.rad = float(fds[9])
            if len(fd) >= 11:
                atom.elem = fds[10]
            else:
                atom.elem = atom.name[0]
        except ValueError:
            print('ATOM line does not conform to pqr standard!\n')
            exit(1)
    else:
        try:
            atom = Atom()
            atom.atid = int(line[6:11])
            atom.name = line[12:16].strip()
            atom.loc = line[16:17]
            atom.resn = line[17:21].strip()
            if line[26:27].isdigit():
                # charmm pdb reserves icode for resi
                atom.resi = int(line[22:27])
            else:
                atom.resi = int(line[22:26])
                atom.icode = line[26]
            atom.chid = line[21:22]
            atom.r = (float(line[30:38]),float(line[38:46]),float(line[46:54]))
            #-------------------------------
            try:
                atom.oc = float(line[54:60])
            except ValueError:
                atom.oc = 1.0
            #-------------------------------
            try:
                atom.bf = float(line[60:66])
            except ValueError:
                atom.bf = 0.0
            #-------------------------------
            atom.sgid = line[72:76]
            if len(line) > 76:
                atom.elem = line[76:78].strip()
            else:
                atom.elem = atom.name[0]
        except ValueError:
            print('ATOM line does not conform to pdb standard!\n')
            exit(1)
    return atom

# ==============================================================================
def writeatom(at):
    fmt = 'ATOM  %5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      ' \
          '%-4s%2s  \n'
    name = (' '+at.name+' '*(3-len(at.name)) if len(at.name)<4 else at.name)
    return fmt%(at.atid, name, at.loc, at.resn, at.chid, at.resi, at.icode,
                at.r[0], at.r[1], at.r[2], at.oc, at.bf, at.sgid, at.elem)

# ==============================================================================
def atdist(at1, at2):
    r1 = at1.r
    r2 = at2.r
    return sqrt((r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2)

# ==============================================================================
def getats(obj):
    """
    Generate an atom list from Pdb, Mol, Segment, or Residue.
    If obj is an atom list, the obj itself will be returned.
    """
    if isinstance(obj, (Segment,Mol,Pdb)):
        return obj.getats()
    elif isinstance(obj, Residue):
        return obj.atoms
    elif isinstance(obj, (tuple,list)):
        if len(obj)==0 or isinstance(obj[0], Atom):
            return obj
        elif isinstance(obj[0], Residue):
            return reduce(add, map(lambda x:x.atoms,obj))
        else:
            return None
    else:
        return None

# ==============================================================================
def getat(obj, resi, name):
    atoms = getats(obj)
    ats = list(filter(lambda x: x.resi==resi and x.name==name, atoms))
    if len(ats) == 1:
        return ats[0]
    elif len(ats) == 0:
        return None
    else:
        print('ERROR (getat): more than 1 atom are found!')
        exit(1)

# ==============================================================================
def sortat(obj):
    """
    Sort atoms in each residue according to atid
    obj: Pdb, Mol, Segment, or res list
    """
    reses = getreses(obj)
    for res in reses:
        res.atoms.sort(key=lambda at: at.atid)

# ==============================================================================
def getreses(obj):
    """
    Generate an res list from Pdb, Mol, Segment.
    If obj is an res list, the obj itself will be returned.
    """
    if isinstance(obj, (Mol,Pdb)):
        return obj.getreses()
    elif isinstance(obj, Segment):
        return obj.reses
    elif isinstance(obj, (tuple,list)):
        return obj
    else:
        return None

# ==============================================================================
def getres(obj, resi, icode=None):
    reses = getreses(obj)
    if icode == None:
        reses = list(filter(lambda x: x.resi==resi, reses))
    elif isinstance(icode, str):
        reses = list(filter(lambda x: x.resi==resi and x.icode==icode, reses))
    else:
        print('ERROR (getres): please input proper icode value')
        exit(1)
    if len(reses) == 1:
        return reses[0]
    elif len(reses) == 0:
        return None
    else:
        print('ERROR (getres): more than 1 residues are found!')
        exit(1)

