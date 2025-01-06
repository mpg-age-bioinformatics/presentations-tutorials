import warnings
warnings.filterwarnings("ignore")
from scipy.stats import pearsonr 
from scipy.stats import spearmanr
from datetime import datetime
from time import process_time 
import multiprocessing as mp
import pandas as pd
import numpy as np
import itertools
import sys
import os
import time

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string, ":: started")

# Start the stopwatch / counter  
t1_start = process_time()  

gct="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"
output="/nexus/posix0/MAGE-flaski/service/posit/home/jboucas/gtexcorr/"
if not os.path.isdir( output ):
    os.makedirs( output )

n_processors=32

def scan_linepos(path):
    """return a list of seek offsets of the beginning of each line"""
    linepos = []
    offset = 0
    with open(path) as inf:     
        # WARNING: CPython 2.7 file.tell() is not accurate on file.next()
        for line in inf:
            linepos.append(offset)
            offset += len(line)
    return linepos

def sample_lines( linepos_input, p, path=gct):
    """return nsamp lines from path where line offsets are in linepos"""
    
    target=5000

    linepos_input=[linepos_input[i:i + target] for i in range(0, len(linepos_input), target)]
    
    # resulst_=[]
    c=0
    for linepos in linepos_input:
        c=c+target

        linepos.sort()  # this may make file reads more efficient

        filename=f"{output}{linepos[0][0]}.{linepos[0][1]}.tsv"

        if os.path.isfile(filename) :
            return ""

        results = []

        linepos=pd.DataFrame(linepos)

        for offset0 in list(set(linepos[0].tolist() )):

            with open(path) as inf:
                inf.seek(offset0)
                l0=inf.readline().split("\n")[0]

            def _corr(offset1,l0=l0):

                with open(path) as inf:

                    inf.seek(offset1)
                    l1=inf.readline().split("\n")[0]

                    l=l0.split("\t")
                    l_=l1.split("\t")
                    gid=l[0]
                    gid_=l_[0]
                    l=l[2:]
                    l_=l_[2:]

                    tmp=pd.DataFrame( {0:l,1:l_} )
                    tmp=tmp[ ( tmp[0]!="" ) & (tmp[1]!="") ]
                    #print(tmp)
                    tmp=tmp.astype(float)
                    tmp=tmp[ ( tmp[0]>0 ) & (tmp[1]>0) ]

                    n=len(tmp)

                    if n > 2 :

                        tmp=np.log10(tmp)

                        l=tmp[0].tolist()
                        l_=tmp[1].tolist()

                        pearson_stat, pearson_p=pearsonr( l, l_ )
                        spearman_corr, spearman_p=spearmanr( l, l_ )

                        res=[ gid, gid_, n, pearson_stat, pearson_p, spearman_corr, spearman_p]

                    else:
                        res=[ gid, gid_, n, None, None, None, None]

                    res=[ str(s) for s in res ] 
                    res="\t".join(res)

                return res

            linepos_=linepos[linepos[0]==offset0]
            linepos_[2]=linepos_[1].apply(lambda x: _corr(x) )

            r0="\n".join( linepos_[2].tolist() )

            results.append(r0)
        
        results="\n".join(results)
        
        # results_.append(results)

        with open(filename, "w") as f:
            f.write(results)
        
        print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", "part:", p,",", c, "done" )
        sys.stdout.flush()
        
        
    # results_="\n".join(results_)
                        
    return p

linepos = scan_linepos(gct) # the scan only need be done once
linepos = linepos[3:] # remove the header lines
combinations = list(itertools.combinations(linepos, 2))
print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", len(combinations), "combinations" )
sys.stdout.flush()

# dev/demo only do the first 40 combinations
# combinations=combinations[:40]

# we want to have chuncks of len 500 each
# target=5000

# lines_to_process=[combinations[i:i + target] for i in range(0, len(combinations), target)]

size_blocks=int(len(combinations)/n_processors)

lines_to_process=[combinations[i:i + size_blocks] for i in range(0, len(combinations), size_blocks)]

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", len(lines_to_process), "parts" )
sys.stdout.flush()

results = []

i=0      

# for lines_to_process_ in lines_to_process :
    
    # i=i + len(lines_to_process_)*target

print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", "starting pool" )
sys.stdout.flush()
    
pool = mp.Pool(n_processors)
for d in lines_to_process:
    output = pool.apply_async(sample_lines, [d,i])
    i=i+1
    results.append(output)
pool.close() # no more tasks
time.sleep(10)
    # print(datetime.now().strftime("%d/%m/%Y %H:%M:%S"), "::", i, "done" )
    # sys.stdout.flush()
    # pool.join()

results=[ s.get() for s in results ] 
results="\n".join(results)
# print(f"Done\n:{results}"

# results=results.split("\n")
# results=[ s.split("\t") for s in results ]
# results=pd.DataFrame(results, columns=[ "gid", "gid_", "n", "pearson_stat", "pearson_p", "spearman_corr", "spearman_p"])
# results.to_csv("gtex.corr.all.tsv", index=None, sep="\t")
# results.to_excel("gtex.corr.all.xlsx", index=None)

# Stop the stopwatch / counter 
t1_stop = process_time() 

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string, ":: finished", )
print("Elapsed time:", t1_stop, t1_start) 
sys.exit(0)
