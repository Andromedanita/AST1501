import numpy            as np
import matplotlib.pylab as plt
import os
import sys
from nemo_funcs import *

sys.path.insert(0, '/Users/anita/dyn-modeling-streams-2014/sim')
#import split_snaps
from split_snaps import split_snaps

def funcc(filename, dir):
    
    '''
    Parameters
    ---------------------------------------------------------------------------------
        filename: location of the nemo output file (all snaps)
        
        dir : directory path of where to copy the splitted snaps files
    
    
    Return
    ---------------------------------------------------------------------------------
        Generate a concatenated version of action-angle-freuqnecy of the final
        snapshot of the nemo simulation
        
    '''
    
    
    # splitting snaps
    print "Splitting snaps started..."
    split_snaps(snapfile=filename)
    print "Done splitting snaps..."
    
    # base name of the outputs
    basefile = filename.split('.')[0]
    filename_dir = filename.split('gd')[0]
    
    # copying the splitted snaps into the specified directory
    print "Copying files..."
    os.system('cp ' +  filename_dir + basefile + '* ' + dir)
    print "Done copying files..."
    
    # reading final snapshot values
    print "Reading the input..."
    mass, pos, vel = nemo_read_output(basefile + '_00041.dat')

    # computing action-angle-frequency of the tail for the final snapshot
    print "Stated 1st set of action-angle-frequencies..."
    val_tail1 = nemo_coord_convert(pos,vel, 0.9, 1.20, True, 8., 220., 0, 5000)
    print "Done first set of action-angle-frequencies..."
    
    print "Stated 2nd set of action-angle-frequencies..."
    val_tail2 = nemo_coord_convert(pos,vel, 0.9, 1.20, True, 8., 220., 5000, 10000)
    print "Done first set of action-angle-frequencies..."


    # saving the action-angle for the tail for both sets of output
    print "Saving the output..."
    np.savetxt("tails_5000.txt" , val_tail1, fmt="%.9e")
    np.savetxt("tails_10000.txt", val_tail2, fmt="%.9e")

    # concatenating the two files fo action-angle
    print "Concatenating output..."
    val_final = np.concatenate((val_tail1, val_tail2), axis=0)

    W0 = float(filename[10:13])
    # saving the concatenated output
    print "Saving final output..."
    np.savetxt("concatenated_{0}.txt".format(W0), val_final, fmt="%.9e")





