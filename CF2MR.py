import re
import numpy as np
periodic_table = {'H': 1.008,
 'He': 4.0026,
 'Li': 6.94,
 'Be': 9.012182,
 'B': 10.81,
 'C': 12.011,
 'N': 14.007,
 'O': 15.999,
 'F': 18.9984032,
 'Ne': 20.1797,
 'Na': 22.98976928,
 'Mg': 24.305,
 'Al': 26.9815386,
 'Si': 28.0855,
 'P': 30.973762,
 'S': 32.06,
 'Cl': 35.45,
 'Ar': 39.948,
 'K': 39.0983,
 'Ca': 40.078,
 'Sc': 44.955912,
 'Ti': 47.867,
 'V': 50.9415,
 'Cr': 51.9961,
 'Mn': 54.938045,
 'Fe': 55.845,
 'Co': 58.933195,
 'Ni': 58.6934,
 'Cu': 63.546,
 'Zn': 65.38,
 'Ga': 69.723,
 'Ge': 72.63,
 'As': 74.9216,
 'Se': 78.96,
 'Br': 79.904,
 'Kr': 83.798,
 'Rb': 85.4678,
 'Sr': 87.62,
 'Y': 88.90585,
 'Zr': 91.224,
 'Nb': 92.90638,
 'Mo': 95.96,
 'Tc': 97.90721,
 'Ru': 101.07,
 'Rh': 102.9055,
 'Pd': 106.42,
 'Ag': 107.8682,
 'Cd': 112.411,
 'In': 114.818,
 'Sn': 118.71,
 'Sb': 121.76,
 'Te': 127.6,
 'I': 126.90447,
 'Xe': 131.293,
 'Cs': 132.9054519,
 'Ba': 137.327,
 'La': 138.90547,
 'Ce': 140.116,
 'Pr': 140.90765,
 'Nd': 144.242,
 'Pm': 144.91276,
 'Sm': 150.36,
 'Eu': 151.964,
 'Gd': 157.25,
 'Tb': 158.92535,
 'Dy': 162.5,
 'Ho': 164.93032,
 'Er': 167.259,
 'Tm': 168.93421,
 'Yb': 173.054,
 'Lu': 174.9668,
 'Hf': 178.49,
 'Ta': 180.94788,
 'W': 183.84,
 'Re': 186.207,
 'Os': 190.23,
 'Ir': 192.217,
 'Pt': 195.084,
 'Au': 196.966569,
 'Hg': 200.59,
 'Tl': 204.383,
 'Pb': 207.2,
 'Bi': 208.9804,
 'Po': 208.98243,
 'At': 209.98715,
 'Rn': 222.01758,
 'Fr': 223.01973,
 'Ra': 226.02541,
 'Ac': 227.02775,
 'Th': 232.03806,
 'Pa': 231.03588,
 'U': 238.02891,
 'Np': 237.04817,
 'Pu': 244.0642,
 'Am': 243.06138,
 'Cm': 247.07035,
 'Bk': 247.07031,
 'Cf': 251.07959,
 'Es': 252.083,
 'Fm': 257.09511,
 'Md': 258.09843,
 'No': 259.101,
 'Lr': 262.11,
 'Rf': 267.122,
 'Db': 268.126,
 'Sg': 271.134,
 'Bh': 274.144,
 'Hs': 277.152,
 'Mt': 278.156,
 'Ds': 281.165,
 'Rg': 282.169,
 'Cn': 285.177,
 'Nh': 286.183,
 'Fl': 289.191,
 'Mc': 290.196,
 'Lv': 293.205,
 'Ts': 294.211}
def separate_chemical_formula(cf):
    m = re.findall('[A-Z][^A-Z]*', cf)
    for i in range(len(m)):
        if not m[i][-1].isdigit():
            m[i] += '1'
    m = ''.join(m)
    el = re.findall('[A-Z][a-z]?', m)
    ar = re.findall('\d+\.?\d*', m)
    ar = [float(x) for x in ar]
    return el, ar
def normal(el0, ar0, el1, ar1):
    ar = [0] * len(ar0)
    for i in range(len(el0)):
        for j in range(len(el1)):
            if el1[j] == el0[i]:
                ar[i] = ar1[j]
    return ar
def solve(list1, list2):
    A = np.array(list1)
    A = A.T
    b = np.array(list2)
    x = np.linalg.lstsq(A, b, rcond=None)
    return x[0]
def prepare(target, part):
    ta0, ta1 = separate_chemical_formula(target)
    pa0 = [[0]*len(ta0)]*len(part)
    pa1 = [ta1]*len(part)
    for i in range(len(part)):
        pa0[i], pa1[i] = separate_chemical_formula(part[i])
    for j in range(len(part)):
        pa1[j] = normal(ta0, ta1, pa0[j], pa1[j])
    return pa1, ta1
def molmass(cf):
    el, ar = separate_chemical_formula(cf)
    molmass = 0
    for i in range(len(el)):
        molmass += periodic_table[el[i]] * ar[i]
    return molmass
target = input('整点儿____\n')
part = input('用____整儿(除字母数字小数点外均可用作分隔符)\n')
part = re.findall('[A-Z][A-Za-z\d\.]*', part)
totalmass = input('整儿____克(g)\n')
totalmass = float(totalmass)
A, b = prepare(target, part)
ratio = solve(A, b)
massratio = [0]*len(part)
partmass = [0]*len(part)
for i in range(len(part)):
    massratio[i] = ratio[i]*molmass(part[i])
    partmass[i] = massratio[i]*totalmass/molmass(target)
ratioprint = [0]*len(part)
for i in range(len(part)):
    ratioprint[i] = str(ratio[i]) + ' ' + part[i]
print(" + ".join(str(i) for i in ratioprint) + ' = ' + target)
massprint = [0]*len(part)
for i in range(len(part)):
    massprint[i] = str(partmass[i])+'g '+ part[i]
print(" + ".join(str(i) for i in massprint) + ' = ' + str(totalmass) + 'g ' + target)
print('质量比为')
mr = [0]*len(part)
for i in range(len(part)):
    mr[i] = [0]*len(part)
    for j in range(len(part)):
        mr[i][j] = massratio[j]/massratio[i]
    print(" : ".join(str(i) for i in part) + ' = ' + " : ".join(str(i) for i in mr[i]))

import msvcrt
print('按任意键结束')
msvcrt.getwche()
