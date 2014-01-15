#this will convert a modflow-usg module into a modflow-2015 module that has
#a derived type with all the variables in it.  note that all variables
#must have the POINTER attribute.  Also, all parameters must be written as
# integer, parameter :: ip=100


import sys

infile = sys.argv[1]
modulename = sys.argv[2]
moduletype = sys.argv[3]
outfile = sys.argv[4]
indent = '  '

print 'infile: ', infile
modulename = modulename.upper()
print 'modulename: ', modulename
print 'outfile: ', outfile

inmodule = False
linelist = []
descripter_var_list = []
f = open(infile, 'r')
for line in f:

    line = line.upper()

    #skip comments
    if line[0] != ' ':
        continue

    if 'MODULE' in line and modulename in line:
        inmodule = True

    if 'END ' in line and modulename in line:
        inmodule = False
        break

    #skip line if not inmodule
    if not inmodule:
        continue

    if line[5] != ' ':
        lastline = linelist.pop()
        line = lastline + line[6:].strip()
        descripter_var_list.pop()
        

    #clean up the line
    line = line.strip()
    
    #create descripter and var
    splitline = line.split('::')
    if len(splitline) == 2:
        descripters = splitline[0].strip().split(',')
        variables = splitline[1].strip().split(',')
        descripter_var_list.append([descripters, variables])

    #save the line    
    linelist.append(line)

f.close()

f = open(outfile, 'w')
line = 'MODULE ' + modulename + '\n'
f.write(line)
line = indent + 'IMPLICIT NONE' + '\n'
f.write(line)


#write the module data
thisindent = indent
for descripters, variables in descripter_var_list:
    descripter = ''
    for d in descripters:
        descripter += d.strip() + ','
    descripter = descripter[0:-1]
    for var in variables:
        newline = thisindent + descripter + ' :: ' + var + '\n'
        f.write(newline)

#write the derived type
line = indent + 'TYPE :: ' + moduletype + '\n'
f.write(line)
thisindent += indent
for descripters, variables in descripter_var_list:
    if 'PARAMETER' in descripters:
        continue
    descripter = ''
    for d in descripters:
        if 'save' in d.lower():
            continue
        descripter += d.strip() + ','
    descripter = descripter[0:-1]
    for var in variables:
        newline = thisindent + descripter + ' :: ' + var + '\n'
        f.write(newline)

#write contained procedures
line = indent + 'CONTAINS' + '\n'
f.write(line)
line = thisindent + 'PROCEDURE :: DESTROY' + '\n'
f.write(line)
line = thisindent + 'PROCEDURE :: PNTSAV' + '\n'
f.write(line)
line = thisindent + 'PROCEDURE :: PNTSET' + '\n'
f.write(line)
thisindent = indent
line = thisindent + 'END TYPE ' + moduletype + '\n'
f.write(line)

#write module contains
#write contained procedures
line = indent + 'CONTAINS' + '\n'
f.write(line)

#write the destroy method
line = indent + 'SUBROUTINE DESTROY(THIS)' + '\n'
f.write(line)
thisindent = 2 * indent
line = thisindent + 'IMPLICIT NONE' + '\n'
f.write(line)
line = thisindent + 'CLASS(' + moduletype + '),INTENT(INOUT) :: THIS' + '\n'
f.write(line)
for descripters, variables in descripter_var_list:
    if 'PARAMETER' in descripters:
        continue
    for var in variables:
        line = thisindent + 'DEALLOCATE(THIS%' + var + ')' + '\n'
        f.write(line)
line = indent + 'END SUBROUTINE DESTROY' + '\n'
f.write(line)

#write the pntsav method
line = indent + 'SUBROUTINE PNTSAV(THIS)' + '\n'
f.write(line)
thisindent = 2 * indent
line = thisindent + 'IMPLICIT NONE' + '\n'
f.write(line)
line = thisindent + 'CLASS(' + moduletype + '),INTENT(INOUT) :: THIS' + '\n'
f.write(line)
for descripters, variables in descripter_var_list:
    if 'PARAMETER' in descripters:
        continue
    for var in variables:
        line = thisindent + 'THIS%' + var + '=>' + var + '\n'
        f.write(line)
line = indent + 'END SUBROUTINE PNTSAV' + '\n'
f.write(line)

#write the pntset method
line = indent + 'SUBROUTINE PNTSET(THIS)' + '\n'
f.write(line)
thisindent = 2 * indent
line = thisindent + 'IMPLICIT NONE' + '\n'
f.write(line)
line = thisindent + 'CLASS(' + moduletype + '),INTENT(INOUT) :: THIS' + '\n'
f.write(line)
for descripters, variables in descripter_var_list:
    if 'PARAMETER' in descripters:
        continue
    for var in variables:
        line = thisindent + var + '=>' +'THIS%' +  var + '\n'
        f.write(line)
line = indent + 'END SUBROUTINE PNTSET' + '\n'
f.write(line)

line = 'END MODULE ' + modulename + '\n'
f.write(line)
