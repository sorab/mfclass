#packages to convert
#bas -- done
#dis -- done
#oc -- done
#bcf/lpf -- done
#wel -- done
#ghb -- done

#sub
#lak
#gag
#sfr

#this script prepares files for converting these
#chd7
#drn7
#evt8
#fhb7
#hfb7
#rch8
#riv7
#str7

from usgmod2mf2015mod import cvtmodule

usgpth = '/Users/langevin/projects/mfusg/mfusg.1_1/src'
outpth = 'createdfiles'
mfusg = usgpth + '/mfusg.f'
cvtlist = []
packlist = [('chd',7), ('drn',7), ('evt',8), ('fhb',7), 
            ('hfb',7), ('rch',8), ('riv',7), ('str',7)]
for p,i in packlist:
    plist = []
    plist.append(usgpth + '/gwf2' + p + str(i) + 'u1.f')
    plist.append('GWF' + p.upper() + 'MODULE')
    plist.append('GWF' + p.upper() + 'TYPE')
    plist.append(outpth + '/gwf3_' + p + 'module.f90')
    cvtlist.append(plist)

#cvtlist = [['gwf2chd7u1', 'GWFCHDMODULE', 'GWFCHDTYPE', 'gwf_chdmodule.f90'],
#           ...
#            ]

#make the module files
for infile, modulename, moduletype, outfile in cvtlist:
    cvtmodule(infile, modulename, moduletype, outfile)

fname = outpth + '/gwf3.txt'
f = open(fname, 'w')

f.write('top of gwf module' + '\n')
for p,i in packlist:
    s = '    type(gwf' + p + 'type) :: gwf' + p + 'dat'
    f.write(s + '\n')
f.write('\n')
    
fusg = open(mfusg,'r')
modelactions = ['ar', 'st', 'rp', 'ad', 'fm', 'bd']
for ma in modelactions:
    f.write('gwf3' + ma + '\n')
    for p,i in packlist:
        s = 'GWF2' + p.upper() + str(i) +'U1' + ma.upper()
        fusg.seek(0)
        for line in fusg:
            if s in line:
                f.write(line.strip() + '\n')
                break
    f.write('\n')
fusg.close()
f.close()
        
        
        