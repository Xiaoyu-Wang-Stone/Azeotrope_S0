import os, sys
import MDAnalysis as mda
import numpy as np 
from math import pi
import time

def read_trj(file_path):
    #os.chdir(file_path)
    u    = mda.Universe('eql.gro', 'eql.trr')
    cell = u.trajectory.ts.dimensions[:3]  # Ångström
    # scaling atomic coordinates 
    sq_time = []
    oxygen_met = u.select_atoms("resname MET and name O*")
    oxygen_water = u.select_atoms("resname SOL and name OW")
    for ts in u.trajectory:
        q1 = oxygen_met.positions
        q2 = oxygen_water.positions
        q_all = np.concatenate((q1, q2), axis=0) 
        sq = np.divide(q_all, cell)
        sq_time.append(sq)
    #set up the name list 
    natoms = oxygen_met.n_atoms + oxygen_water.n_atoms
    names = np.zeros(natoms,'U3')
    names[:oxygen_met.n_atoms] = "Met"
    names[oxygen_met.n_atoms:] = "SOL"
    return [cell, names, sq_time]


def FT_density(q, kgrid):
    # This is the un-normalized FT for density fluctuations
    ng = len(kgrid)
    ak = np.zeros(ng,dtype=complex)
    for n,k in enumerate(kgrid):
        ak[n] = np.sum(np.exp(-1j*(q[:,0]*k[0]+q[:,1]*k[1]+q[:,2]*k[2])))
    return ak

def Sk(names, q, kgrid, e_A, e_B):
    # This is the un-normalized FT for the density fluctuations
    q_A = np.asarray([ q_now for i, q_now in enumerate(q) if names[i] in e_A ])
    n_A = len(q_A)
    print("Number of element A: ", n_A)
    if n_A > 0:
        FTrho_A = FT_density(q_A, kgrid)
    else:
        FTrho_A = np.empty(len(kgrid))
        FTrho_A[:] = np.NaN
    if e_A != e_B:
        q_B = np.asarray([ q_now for i,q_now in enumerate(q) if names[i] in e_B ])
        n_B = len(q_B)
        print("Number of element B: ", n_B)
        if n_B > 0:
            FTrho_B = FT_density(q_B, kgrid)
        else:
            FTrho_B = np.empty(len(kgrid))
            FTrho_B[:] = np.NaN
    else:
        FTrho_B = FTrho_A

    return np.multiply(FTrho_A, np.conjugate(FTrho_A))/n_A, \
                   np.multiply(FTrho_A, np.conjugate(FTrho_B))/(n_A*n_B)**0.5, \
                   np.multiply(FTrho_B, np.conjugate(FTrho_B))/n_B



def main(sprefix="Sk", straj="read", sbins=8):
    # the input file
    #os.chdir(straj)
    # number of k grids
    #sbins = 8
    bins = int(sbins)
    # get total number of bins and initialize the grid
    u    = mda.Universe('eql.gro', 'eql.trr')
    print("Use number of bins:", bins)
    
    #sprefix = "Sk"
    # Outputs
    ofile_AA = open(sprefix+'-II-real.dat',"ab")
    ofile_AB = open(sprefix+'-IW-real.dat',"ab")
    ofile_BB = open(sprefix+'-WW-real.dat',"ab")
    
    [cell, names, sq_data] = read_trj(straj)
    
    nframe = 0
    for i in range(len(u.trajectory)):
        start_time = time.time()
        #[cell, names, sq_data] = read_trj(straj)
        nframe += 1
        print("Frame No:", nframe)
    
        if (nframe == 1):
            # normalization
            volume = np.prod(cell[:])
            kgrid = np.zeros((bins*bins*bins,3),float)
            kgridbase = np.zeros((bins*bins*bins,3),float)
            # initialize k grid
            [ dkx, dky, dkz ] = [ 1./cell[0], 1./cell[1], 1./cell[2] ]
            n=0
            for i in range(bins):
                for j in range(bins):
                    for k in range(bins):
                        if i+j+k == 0: pass
                        # initialize k grid
                        kgridbase[n,:] = (2.*pi)*np.array([i, j, k])
                        kgrid[n,:] = [ dkx*i, dky*j, dkz*k ]
                        n+=1
            np.savetxt(sprefix+'-kgrid.dat',kgrid)
        #print("--- %s seconds after read frame ---" % (time.time() - start_time))
        # FT analysis of density fluctuations
        sk_AA, sk_AB, sk_BB = Sk(names, sq_data[i], kgridbase, ['Met'], ['SOL'])
        print("--- %s seconds after FFT density ---" % (time.time() - start_time))
        # Outputs
        np.savetxt(ofile_AA,sk_AA[None].real, fmt='%4.4e', delimiter=' ',header="Frame No: "+str(nframe))
        np.savetxt(ofile_AB,sk_AB[None].real, fmt='%4.4e', delimiter=' ',header="Frame No: "+str(nframe))
        np.savetxt(ofile_BB,sk_BB[None].real, fmt='%4.4e', delimiter=' ',header="Frame No: "+str(nframe))
    
    print("A total of data points ", nframe)
    sys.exit()

if __name__ == '__main__':
    file_loc = sys.argv[0]
    main("Sk",file_loc, 8)    

