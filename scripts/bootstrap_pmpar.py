"""
you can streamline it with it code. the functionality can be tested in ipython.
1. import bootstrap_pmpar.py
2. run bootstrap_pmpar.bootstrap_pmpar()\
    according to its usage help (bootstrap_pmpar will call the two other functions).
"""
import os, sys
import numpy as np
from astropy.table import vstack, Table
def dms2deg(string): #-- convert dd:mm:ss.ssss to dd.dddd
    a = string.split(':')
    b = []
    for item in a:
        b.append(float(item))
    if b[0] < 0:
        sign = -1
    else:
        sign = 1
    b[0] = abs(b[0])
    i=len(b)
    degree = b[i-1]
    while i != 1:
        degree /= 60
        i -= 1
        degree += b[i-1]
    degree *= sign
    return degree
def readpmparout(pmparout):
    rchsq = 0
    lines = open(pmparout).readlines()
    for line in lines:
        if 'epoch' in line:
            epoch = line.split('=')[1].strip()
        if 'Reduced' in line:
            rchsq = float(line.split('=')[1].strip())
        for estimate in ['l', 'b']:
            if (estimate in line) and 'degrees' in line:
                exec("%s = %s" % (estimate, line.split('=')[-1].strip().split(' ')[0]))
        for estimate in ['RA', 'Dec  ', 'mu_a', 'mu_d', 'pi']:
            if estimate in line:
                print line.split('=')[-1].split('+')[0].strip()
                exec("%s = line.split('=')[-1].split('+')[0].strip()" % estimate.strip())
                if estimate in ['RA', 'Dec  ']:
                    exec("%s = dms2deg(%s)" % (estimate.strip(), estimate.strip()))
    for line in lines:
        for estimate in ['mu_a', 'mu_d', 'pi']:
            if estimate in line:
                error = line.split('+-')[1].strip().split(' ')[0]
                exec("error_%s = %s" % (estimate, error))
                exec("%s = float(%s)" % (estimate, estimate))
    return RA, Dec, epoch, pi, mu_a, mu_d, error_pi, error_mu_a, error_mu_d, l, b, rchsq
def pulsitions2paras(pulsitions):
    pmparout = './.pmpar.out'
    os.system("pmpar %s > %s" % (pulsitions, pmparout))
    [RA0, Dec0, junk1, PI0, mu_a0, mu_d0, junk2, junk3, junk4, junk5, junk6, rchsq] = readpmparout(pmparout)
    D0 = 1/PI0
    return D0, PI0, mu_a0, mu_d0, RA0, Dec0, rchsq
def bootstrap_pmpar(pmparinfile, bootstrapruns, priors='', overwrite_table=False):
    """
    1. pmparinfile is the input pmpar.in file;
    2. after running this module, a file called .five_parameters.dat would be generated, which stores the parallaxes, proper motions,\
        and reference positions results from pmpar, and can be written anytime;
    3. if you want to build your bootstrap sample with several runs, set overwrite_table to True.
    """
    pulsitions = './.pmpar.in.bootstrap' 
    if not os.path.exists(pmparinfile):
        print("%s does not exists; aborting\n" % pmparinfile)
        sys.exit()
    positions = []
    lines = open(pmparinfile).readlines()
    for line in lines:
        if 'epoch' in line:
            epochline = line
        if line.count(':') == 4 and (not line.strip().startswith('#')):
            positions.append(line)
    nepoch = len(positions)
    PIs = np.array([])
    mu_as = np.array([])
    mu_ds = np.array([])
    RAs = np.array([])
    Decs = np.array([])
    count = 0
    while count < bootstrapruns:
        random_indices = np.random.randint(0, nepoch, nepoch)
        if len(np.unique(random_indices)) < 3: #use at least 3 different positions for astrometric fit
            continue
        fileWrite = open(pulsitions, 'w')
        fileWrite.write(epochline)
        fileWrite.write(priors)
        for i in random_indices:
            fileWrite.write(positions[i])
        fileWrite.close()
        [D, PI, mu_a, mu_d, RA, Dec, rchsq] = pulsitions2paras(pulsitions) 
        PIs = np.append(PIs, PI)
        mu_as = np.append(mu_as, mu_a)
        mu_ds = np.append(mu_ds, mu_d)
        RAs = np.append(RAs, RA)
        Decs = np.append(Decs, Dec)
        print("\x1B[1A\x1B[Kprogress:{0}%".format(round((count + 1) * 1000 / bootstrapruns)/10) + " \r")
        count += 1
    t = Table([PIs, mu_as, mu_ds, RAs, Decs], names=['PI', 'mu_a', 'mu_d', 'RA', 'Dec']) 
    print t
    bootstrapped_five_parameters_table = './.five_parameters.dat'
    if not overwrite_table:
        if os.path.exists(bootstrapped_five_parameters_table):
            t0 = Table.read(bootstrapped_five_parameters_table, format='ascii')
            t = vstack([t0, t])
    t.write(bootstrapped_five_parameters_table, format='ascii', overwrite=True)
