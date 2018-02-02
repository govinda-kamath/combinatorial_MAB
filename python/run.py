import multiprocessing as mp
import numpy as np
import subprocess
import sys


def callKNN(command):
    subprocess.call(command.split())
    # print command

numThreads = int(sys.argv[1])
maxPoints = float(sys.argv[2])

chunkSize = np.ceil(maxPoints/float(numThreads))

cmdVec = []
for i in np.arange(int(maxPoints), step=int(chunkSize)):
    start = i
    end = i + int(chunkSize)

    cmd = './build/test_code  nominal.ini '+ str(start) + " " + str(end)
    cmdVec.append(cmd)

pool=mp.Pool(processes=numThreads)
pool.map(callKNN,cmdVec)


