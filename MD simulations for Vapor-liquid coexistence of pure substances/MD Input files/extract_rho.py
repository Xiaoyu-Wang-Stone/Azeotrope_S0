#!/bin/bash
import os, shutil, subprocess
import numpy as np

press = []
temper = [330,340,350,360, 380, 390, 400,410 ]

for temp in temper:
    press = []
    for i in range(1,11):
        file_path = f'/nfs/scistore14/chenggrp/xwang/eth_kbff_NVT/{temp}K/prd{i}'
        os.chdir(file_path)

        #command = "gmx energy -f eql.edr -b 200 <<< \"32 0\" | grep -i  \"Pres-ZZ\" | awk \'{print $2 \"\t\" $3}\' > pressZ.txt"
        #subprocess.run(command, shell=True, executable='/bin/bash')
        data = np.loadtxt(f"/nfs/scistore14/chenggrp/xwang/eth_kbff_NVT/{temp}K/prd{i}/pressZ.txt", delimiter='\t')
        press.append([i, float(data[0]), float(data[1])])
    press_data = np.array(press)
    print(f"{temp} {np.mean(press_data[:,1]):2f}  {np.std(press_data[:,1]):2f}")

press_data = np.array(press)

print(press_data)
print(f"{np.mean(press_data[:,1]):2f}")
print(f"{np.std(press_data[:,1]):2f}")
print(f"{np.mean(press_data[:,2]):2f}")
