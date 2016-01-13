import numpy            as np
import matplotlib.pylab as plt
import os

sys.path.insert(0, '/Users/anita/dyn-modeling-streams-2014/sim')
import split_snaps

def (filename, dir):
    
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
    split_snaps(snapfile=filename)
    
    # base name of the outputs
    basefile = filename.split('.')[0]
    filename_dir = filename.split('gd')[0]
    
    # copying the splitted snaps into the specified directory
    os.system('cp ' +  filename_dir + basefile + '* ' + dir)
    
    # reading final snapshot values
    mass, pos, vel = nemo_read_output(basefile + '_00041.dat')

    # computing action-angle-frequency of the tail for the final snapshot
    val_tail1 = nemo_coord_convert(pos,vel, 0.9, 1.20, True, 8., 220., 0, 5000)
    val_tail2 = nemo_coord_convert(pos,vel, 0.9, 1.20, True, 8., 220., 5000, 10000)

    # saving the action-angle for the tail for both sets of output
    np.savetxt("tails_5000.txt" , val_tail1, fmt="%.9e")
    np.savetxt("tails_10000.txt", val_tail2, fmt="%.9e")

    # concatenating the two files fo action-angle
    val_final = np.concatenate((val_tail1, val_tail2), axis=0)

    W0 = float(filename[10:13])
    # saving the concatenated output
    np.savetxt("concatenated_{0}.txt".format(W0), val_final, fmt="%.9e")





