
__doc__="""
Count ions
"""

import sys,os
import random as rand
import numpy as np
from pmx import *
from pmx import library
from pmx import geometry
from pmx import xtc
from pmx import ndx
from pmx.options import *
from pmx.parser import *

mainchain = [ 'N','NH','C','O','CA','OC1','OC2','OX','H','HA','HA1','HA2' ]
bckb = [ 'N', 'C', 'CA', 'NH' ]

def track_permeations( ions, trjTimes, bCyl=False, ionsCyl={} ):
    transitionsUp = {}
    transitionsDown = {}
    transitionsUpTimes = {}
    transitionsDownTimes = {}
    totalUp = 0
    totalDown = 0
    # counting within cylinder (optional)
    CYLtransitionsUp = {}
    CYLtransitionsDown = {}
    CYLtransitionsUpTimes = {}
    CYLtransitionsDownTimes = {}
    CYLtotalUp = 0
    CYLtotalDown = 0

    for key in ions.keys():
        transitionsUp[key] = 0 
        transitionsDown[key] = 0 
        transitionsUpTimes[key] = []
        transitionsDownTimes[key] = []
        # counting within cylinder (optional)
        CYLtransitionsUp[key] = 0 
        CYLtransitionsDown[key] = 0 
        CYLtransitionsUpTimes[key] = []
        CYLtransitionsDownTimes[key] = []

        state0 = 0 # for tracking 1->2->3->4
        state1 = 5 # for tracking 4->3->2->1
        CYLstate0 = 0 # for tracking 2->3 in the cylinder
        CYLstate1 = 0 # for tracking 3->2 in the cylinder

        counter = 0
        for i in ions[key]:
            if i==state0+1: # step forward
                state0+=1
            elif state0==4 and i==1: # periodic step forward
                state0=1
            elif i==state0-1: # step backward
                state0-=1
            elif state0==1 and i==4: # periodic step backward
                state0=0
            # WARNING: jump over two positions
            elif (state0==2 and i==4):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 0
            elif (state0==4 and i==2):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 2
            elif (state0==3 and i==1):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 4
            elif (state0==1 and i==3):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 0
            ######## cylinder #######
            if bCyl==True:
                if state0==2:
                    if ionsCyl[key][counter]==1:
                        CYLstate0 = 1
                    else:
                        CYLstate0 = 0
                if state0==3:
                    if ionsCyl[key][counter]==1 and (CYLstate0==1 or CYLstate0==2):
                        CYLstate0 = 2
                    else:
                        CYLstate0 = 0
            #########################

            if i==state1-1: # step forward
                state1-=1
            elif state1==1 and i==4: # periodic step forward
                state1=4
            elif i==state1+1: # step backward
                state1+=1
            elif state1==4 and i==1: # periodic step backward
                state1=5
            elif (state1==2 and i==4):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 1
            elif (state1==4 and i==2):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 5
            elif (state1==3 and i==1):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 5
            elif (state1==1 and i==3):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 3

            # WARNING: jump over two positions
            elif (state1==2 and i==4) or \
                 (state1==4 and i==2) or \
                 (state1==3 and i==1) or \
                 (state1==1 and i==3):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1=5
            ######## cylinder #######
            if bCyl==True:
                if state1==3:
                    if ionsCyl[key][counter]==1:
                        CYLstate1 = 1
                    else:
                        CYLstate1 = 0
                if state1==2:
                    if ionsCyl[key][counter]==1 and (CYLstate1==1 or CYLstate1==2):
                        CYLstate1 = 2
                    else:
                        CYLstate1 = 0
            #########################
           
            if state0==4: # transition up
                state0 = 0
                transitionsUp[key]+=1
                transitionsUpTimes[key].append(trjTimes[counter])
                totalUp+=1
                if bCyl==True and CYLstate0==2:
                    CYLstate0 = 0
                    CYLtransitionsUp[key]+=1
                    CYLtransitionsUpTimes[key].append(trjTimes[counter])
                    CYLtotalUp+=1    
            if state1==1: # transition down
                state1 = 5
                transitionsDown[key]+=1
                transitionsDownTimes[key].append(trjTimes[counter])
                totalDown+=1
                if bCyl==True and CYLstate1==2:
                    CYLstate1 = 0
                    CYLtransitionsDown[key]+=1
                    CYLtransitionsDownTimes[key].append(trjTimes[counter])
                    CYLtotalDown+=1    

            counter+=1

    return(transitionsUp,transitionsDown,transitionsUpTimes,transitionsDownTimes,totalUp,totalDown,
           CYLtransitionsUp,CYLtransitionsDown,CYLtransitionsUpTimes,CYLtransitionsDownTimes,CYLtotalUp,CYLtotalDown)

def identify_range( z, plow, pmid, phigh ):
    r = 0
    if z < plow:
        r = 1
    elif z>=plow and z<pmid:
        r = 2
    elif z>=pmid and z<phigh:
        r = 3
    elif z>=phigh:
        r = 4
    return(r)

def identify_if_in_cylinder( x, y, cylCenter, cylRad ):
    d = np.sqrt( np.power(x-cylCenter[0],2.0) + np.power(y-cylCenter[1],2.0) )
    if d <= cylRad:
        return(1)
    return(0)

def get_cylinder( crd ):
    # return the [x,y] coordinates of the center point and radius
    center = [0.0,0.0]
    rad = 0.0

    # 1. get the center
    count = 0
    for c in crd:
        center[0] += c[0]
        center[1] += c[1]
        count+=1
    center[0] /= float(count)
    center[1] /= float(count)

    # 2. get the radius (as a maximal distance from the center)
    for c in crd:
        d = np.sqrt( np.power(c[0]-center[0],2.0) + np.power(c[1]-center[1],2.0) )
        if d > rad:
            rad = d

    return(center,rad)

def get_limits( crd, z=2 ):
    pmid = np.mean( crd, axis=0 )[z]
    indhigh = np.where(crd[:,z]>pmid)
    indlow = np.where(crd[:,z]<pmid)
    phigh = np.mean(crd[indhigh,z])
    plow = np.mean(crd[indlow,z])
    return(plow,pmid,phigh)

def select_ndx( fname="index.ndx", message=False ):
    if not os.path.isfile( fname ): return(False)
    ndx_file = ndx.IndexFile( fname )
    names = ndx_file.names
    i = 0
    ndxDict = {}
    for name in names:
        atomNum = len(ndx_file[name].ids)
        sys.stdout.write('%d %s: %d atoms\n' % (i,name,atomNum) )
        ndxDict[i] = ndx_file[name].ids
        i+=1
    sys.stdout.write('\n')

    if message==False:
        sys.stdout.write('Select a group for analysis:\n')
    else:
        sys.stdout.write(message+'\n')

    ndxNum = -1
    while ndxNum==-1:
        ndxNum = input()
        if ndxNum.isdigit()==False:
            sys.stdout.write('Wrong index group number selected (use an integer number)\n')
            ndxNum = -1
            continue
                
        ndxNum = int(ndxNum) 
        if (ndxNum >= i) or (ndxNum < 0):
            sys.stdout.write('Wrong index group number selected\n')
            ndxNum = -1
    sys.stdout.write('Selected group %d\n\n' % ndxNum)

    res = []
    res = np.asarray(ndxDict[ndxNum])-1 # starting from 0
    return(res,ndxNum)
 
    
def main(argv):

   options = [
        Option( "-cyl", "bool", False, "should counting in the cylinder be performed"),
#       Option( "-ncomp", "int", 1, "number of components to use; -1 determines optimal number using cross-validation data; -2 determines optimal number using training data"),
#       Option( "-cv_percent", "float", 0.0, "a number in the interval (0%;100%) will perform cross-validation using this percent of frames for crossval"),
#       Option( "-fit1", "string", "standard", "perform fitting of group 1: standard, bckbaa, none"),
        ]
    
   files = [
       FileOption("-s", "r",["pdb","gro"],"init.pdb", "input structure file"),
       FileOption("-f", "r",["pdb","xtc"],"traj.xtc", "input trajectory file"),
       FileOption("-n", "r",["ndx"],"index.ndx", "input index file"),
       FileOption("-o", "w",["dat"],"out.dat", "output counts"),
       ]

   help_text = ('Count ions',
                'Takes into account all ion crossings from one side of the bilayer to the other.',
                'Considers only one ion species at a time.',
                'To count several ion species run script several times separately.',
                '',
                'Optionally, can consider ions passing through a cylinder.',
                'The cylinder is constructed from the atoms in a separate index group, e.g. protein atoms.',
                'These atoms are used to define height and radius of the cylinder parallel to z-axis.',
                'For the -cyl option to work properly, the system needs to be superimposed accordingly.',)
 

    
   cmdl = Commandline( argv, options = options,
                       fileoptions = files,
                       program_desc = help_text,
                       check_for_existing_files = False )

################################
# scheme for region counting
################################
################################
#
# (bulk) 4 
# PPPP  phigh  PPPPPPPPPPPPPPPPP
#  3
#------  middle  -------------
#  2 
# PPPP  plow   PPPPPPPPPPPPPPPP
# (bulk) 1
#
##################################
##################################
# when counting within the cylinder #
# only regions 2 and 3 need to be within its radius #
#####################################################

   plow = 0.0
   pmid = 0.0
   phigh = 0.0

   #### files ####
   strFile = cmdl['-s']
   trjFile = cmdl['-f']
   ndxFile = cmdl['-n']
   predFile = cmdl['-o']
   bCyl = cmdl['-cyl']

   init = Model(strFile,bPDBTER=True,renumber_atoms=False,renumber_residues=False)
   struct = Model(strFile,bPDBTER=True,renumber_atoms=False,renumber_residues=False)
   trj = xtc.Trajectory(trjFile)

   #### index ####
   ndxPhos = []
   ndxIons = [] 
   ndxCyl = []
   (ndxPhos,ndxPhosNum) = select_ndx(ndxFile,message='Select phosphate group\n')
   (ndxIons,ndxIonsNum) = select_ndx(ndxFile,message='Select ion group\n')
   if bCyl==True:
       (ndxCyl,ndxCylNum) = select_ndx(ndxFile,message='Select protein group for cylinder\n')
   ionNum = len(ndxIons)

   ### initialize ion structures ###
   ions = {}
   ionsCyl = {} # track wheter within a cylinder (optional): 0 - not in cylinder, 1 - in cylinder
   for i in range(0,ionNum):
       ions[ndxIons[i]] = []
       ionsCyl[ndxIons[i]] = []

   ##############################
   #### reading trajectory ####
   frnum = 0
   trjTimes = []
   for frame in trj:                        # go over each frame 
       frame.update( struct )    		    # update coords in model 

       # identify 1-2-3-4 regions and, if needed, cylinder (for each frame)
       if frnum>=0:
           phosCrd = np.array( list(map(lambda i: struct.atoms[i].x, ndxPhos)) )
           plow,pmid,phigh = get_limits( phosCrd )
           if bCyl==True:
               cylCrd = np.array( list(map(lambda i: struct.atoms[i].x, ndxCyl)) )
               cylCenter,cylRad = get_cylinder( cylCrd )
#           print(plow,pmid,phigh)
#           print(cylCenter,cylRad)

       # analyze ions
       for i in range(0,ionNum):
           z = struct.atoms[ndxIons[i]].x[2]
           # identify the range
           r = identify_range( z, plow, pmid, phigh )
           ions[ndxIons[i]].append(r)
           # track if in cylinder
           if bCyl==True:
               x = struct.atoms[ndxIons[i]].x[0]
               y = struct.atoms[ndxIons[i]].x[1]
               c = identify_if_in_cylinder( x, y, cylCenter, cylRad )
               ionsCyl[ndxIons[i]].append(c)
       
       # store trj time
       trjTimes.append(frame.time)

       sys.stdout.write('Step: %d       Time: %f                  Frame: %d\r' % (frame.step,frame.time,frnum))
       sys.stdout.flush()
       frnum+=1
   sys.stdout.write('\n')
   trj.xdr.xdrfile_close(trj.xd) # close


   # track permeations
   transitionsUp,transitionsDown,transitionsUpTimes,transitionsDownTimes,totalUp,totalDown,CYLtransitionsUp,CYLtransitionsDown,CYLtransitionsUpTimes,CYLtransitionsDownTimes,CYLtotalUp,CYLtotalDown = track_permeations( ions, trjTimes, bCyl, ionsCyl )

   ######################################################
   ######################################################
   ###################### output ########################
   ######################################################
   ######################################################
   fp = open(cmdl['-o'],'w')
   # summary
   fp.write('-----------------------------------\n')
   fp.write('-----------------------------------\n')
   fp.write('Total simulation time: {0} ps\n'.format(trjTimes[-1]))
   fp.write('Total permeations up: {0}\n'.format(totalUp))
   fp.write('Total permeations down: {0}\n'.format(totalDown))
   fp.write('-----------------------------------\n')
   currentUp = totalUp*1.602176634/trjTimes[-1]*100000.0 # pA
   currentDown = totalDown*1.602176634/trjTimes[-1]*100000.0 # pA
   fp.write('Current up: {0} pA\n'.format( round(currentUp,5)))
   fp.write('Current down: {0} pA\n'.format( round(currentDown,5)))
   fp.write('-----------------------------------\n')
   fp.write('-----------------------------------\n')

   #########################
   ##### cylinder ##########
   #########################
   if bCyl==True:
       fp.write('\n**************************\n')
       fp.write('******** Cylinder ********\n')
       fp.write('Cylinder permeations up: {0}\n'.format(CYLtotalUp))
       fp.write('Cylinder permeations down: {0}\n'.format(CYLtotalDown))
       fp.write('---------\n')
       currentUp = CYLtotalUp*1.602176634/trjTimes[-1]*100000.0 # pA
       currentDown = CYLtotalDown*1.602176634/trjTimes[-1]*100000.0 # pA
       fp.write('Cylinder current up: {0} pA\n'.format( round(currentUp,5)))
       fp.write('Cylinder current down: {0} pA\n'.format( round(currentDown,5)))
       fp.write('******** Cylinder ********\n')
       fp.write('**************************\n\n')

   fp.write('-----------\n')
   fp.write('--Details--\n')
   fp.write('-----------\n')
   if totalUp>0:
       for key in transitionsUp.keys():
           if transitionsUp[key]>0:
               fp.write('Up: ion {0} had {1} permeations at times (ps):'.format(key,transitionsUp[key]))
               for i in range(0,transitionsUp[key]):
                   fp.write(' {0}'.format(transitionsUpTimes[key][i]))
               fp.write('\n')
   if totalDown>0:
       for key in transitionsDown.keys():
           if transitionsDown[key]>0:
               fp.write('Down: ion {0} had {1} permeations at times (ps):'.format(key,transitionsDown[key]))
               for i in range(0,transitionsDown[key]):
                   fp.write(' {0}'.format(transitionsDownTimes[key][i]))
               fp.write('\n')

   #########################
   ##### cylinder ##########
   #########################
   if bCyl==True:
       fp.write('\n***********\n')
       fp.write('**Details**\n')
       fp.write('***********\n')
       if CYLtotalUp>0:
           for key in CYLtransitionsUp.keys():
               if CYLtransitionsUp[key]>0:
                   fp.write('Cylinder up: ion {0} had {1} permeations at times (ps):'.format(key,CYLtransitionsUp[key]))
                   for i in range(0,CYLtransitionsUp[key]):
                       fp.write(' {0}'.format(CYLtransitionsUpTimes[key][i]))
                   fp.write('\n')
       if CYLtotalDown>0:
           for key in CYLtransitionsDown.keys():
               if CYLtransitionsDown[key]>0:
                   fp.write('Cylinder down: ion {0} had {1} permeations at times (ps):'.format(key,CYLtransitionsDown[key]))
                   for i in range(0,CYLtransitionsDown[key]):
                       fp.write(' {0}'.format(CYLtransitionsDownTimes[key][i]))
                   fp.write('\n')
   fp.close()

#   for i in range(0,ionNum):
#       if ndxIons[i]==30976:
#           for j in range(0,np.shape(ions[ndxIons[i]])[0]):
#               print(j,ndxIons[i],ions[ndxIons[i]][j],ionsCyl[ndxIons[i]][j])
#       print(ndxIons[i],ionsCyl[ndxIons[i]])
   

if __name__=='__main__':
    main(sys.argv)

