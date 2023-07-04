import os
import shutil
from hatchet.utils.Supporting import ensure, log, error
from util_pcf import forward,backward,central,estimate_forward



def select(
    diploid,
    tetraploid,
    rundir="",
    g=0.3,
    limit=0.6
):
    assert len(diploid) == len(tetraploid), error('The number of diploid and tetraploid results must be the same')
    dscores={}
    tscores={}
    clones=list(diploid.keys())
    if len(clones) == 1 or len(clones) == 2:
        dscores=diploid
        tscores=tetraploid
    else:
        for i in clones:
            if i == min(clones):
                dscores[i] = forward(diploid, i, g, limit)
            elif i == max(clones):
                dscores[i] = backward(diploid, i, g, limit)
            else:
                dscores[i] = central(diploid, i, g, limit)


        for i in clones:
            if i == min(clones):
                tscores[i] = estimate_forward(tetraploid, i, g, diploid, limit)
            elif i == max(clones):
                tscores[i] = backward(tetraploid, i, g, limit)
            else:
                tscores[i] = central(tetraploid, i, g, limit)


    dchosen = max(diploid, key=(lambda x: dscores[x]))
    dbout = os.path.join(rundir, 'chosen.diploid.bbc.ucn')
    shutil.copy2("results.diploid.n"+str(dchosen) + '.bbc.ucn.tsv', dbout)

  
    tchosen = max(tetraploid, key=(lambda x: tscores[x]))
    tbout = os.path.join(rundir, 'chosen.tetraploid.bbc.ucn')
    shutil.copy2(os.path.join(rundir,"results.tetraploid.n"+str(tchosen) + '.bbc.ucn.tsv'), tbout)
  

    bbest = os.path.join(rundir, 'best.bbc.ucn')
    
    if tchosen < dchosen:
        shutil.copy2(os.path.join(rundir,"results.diploid.n"+str(tchosen) + '.bbc.ucn.tsv'), bbest)
    
    else:
        shutil.copy2(os.path.join(rundir,"results.diploid.n"+str(dchosen) + '.bbc.ucn.tsv'), bbest)

    return dchosen,tchosen
