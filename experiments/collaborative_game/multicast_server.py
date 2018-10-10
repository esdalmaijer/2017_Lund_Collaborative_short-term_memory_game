from constants import *
from pair import PairServer

import copy
import time
import random
from itertools import combinations    
from multiprocessing import Process
    
# FUNCTION
def start_pair(ip_addresses, nstim, rewtypes, ntrials):
    
    """Starts a Server class for a pair of client numbers.
    """

    pair = PairServer(ip_addresses)
    pair.prepare_block(nstim, rewtypes, ntrials)
    pair.run_block()
    pair.close()
    
    return 0


if __name__ == "__main__":

    print("Running!")
    
    # PREPARATION
    # Manually create a unique set of pairs, as the programmatic way is
    # a fucking hassle.
    if len(CLIENTIPS) == 2:
        all_pairs = [ \
            [[1, 2]] \
            ]
    elif len(CLIENTIPS) == 4:
        all_pairs = [ \
            [[1,2], [3,4]], \
            [[1,3], [2,4]], \
            [[1,4], [2,3]] \
            ]
    elif len(CLIENTIPS) == 6:
        all_pairs = [ \
            [[1,2], [3,4], [5,6]], \
            [[1,3], [2,5], [4,6]], \
            [[1,4], [2,6], [3,5]], \
            [[1,5], [2,4], [3,6]] \
#            [[1,6], [2,3], [4,5]] \
            ]
    elif len(CLIENTIPS) == 8:
        all_pairs = [ \
            [[1,2], [3,4], [5,6], [7,8]], \
            [[1,3], [2,4], [5,7], [6,8]], \
            [[1,4], [2,3], [5,8], [6,7]], \
            [[1,5], [2,6], [3,7], [4,8]] \
#            [[1,6], [2,5], [3,8], [4,7]], \
            ]
    elif len(CLIENTIPS) == 10:
        all_pairs = [ \
            [[1,2], [3,4], [5,6], [7,8], [9,10]], \
            [[1,3], [2,4], [5,7], [6,9], [8,10]], \
            [[1,4], [2,3], [5,8], [6,10], [7,9]], \
            [[1,5], [2,6], [3,7], [4,10], [8,9]] \
            ]

    nblocks = len(all_pairs)
#
#    nblocks = len(CLIENTIPS) - 1
#    
#    # Generate unique pairs of IPs.
#    combs = combinations(CLIENTIPS, 2)
#    all_pairs = list(combs)
#    blocks = []
#    while len(all_pairs) > 0:
#        print("Generating a new block")
#        current_ips = []
#        thisblock = []
#        while len(current_ips) < len(CLIENTIPS):
#            i = 0
#            while (all_pairs[i][0] in current_ips) \
#                or (all_pairs[i][1] in current_ips):
#                i = random.randint(0, len(all_pairs)-1)
#            p = all_pairs.pop(i)
#            print("\tSelected pair %s" % list(p))
#            current_ips.extend(list(p))
#            thisblock.append(copy.deepcopy(p))
#        print("Adding new block to the list!")
#        blocks.append(copy.deepcopy(thisblock))

    # EXPERIMENT
    for blocknr in range(nblocks):
    
        # Create pairs from IPs
        ip_addresses = []
        for i in range(len(all_pairs[blocknr])):
            ip_addresses.append([CLIENTIPS[all_pairs[blocknr][i][0]-1], \
                CLIENTIPS[all_pairs[blocknr][i][1]-1]])
        
        # Run each pair in a parallel Process
        processes = []
        for ip_pair in ip_addresses:
            
            p = Process( \
                target=start_pair, \
                args=(ip_pair, NSTIM, REWTYPES, NTRIALS) \
                )
            p.name = "pair_%s" % ('_'.join(map(str,ip_pair)))
            processes.append(p)
        
        for p in processes:
            print("Starting process %s" % (p.name))
            p.start()
        
        for p in processes:
            if p.is_alive():
                print("Waiting for process %s to finish..." % (p.name))
                p.join()
            else:
                print("Process %s died an unexpected death... RIP." % (p.name))
                p.terminate()
    
        print("All pairs are finished.")
    
        # Wait for keyboard input to allow people to answer the post-run questions.
        if blocknr < nblocks - 1:
            os.system("pause")
