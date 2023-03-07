import os
import shutil
import sys
import subprocess
import parmed as pmd
import argparse
from multiprocessing import Pool
from tqdm import tqdm
import socket

def write_gauss_shell(curr):
    # Making Shell Script Portable -- only usable for Fourier RN
    gauss='''#!/bin/bash
    jobn=$1
    outn=${jobn//".com"/}
    export GAUSS_SCRDIR=/mnt/data1/g09scr
    export GAUSS_EXEDIR=/opt/g09E/
    mkdir -p $GAUSS_SCRDIR
    /opt/g09E/g09 <$jobn > $outn.log'''
    with open("{0}/gauss.sh".format(curr),"w") as f:
        f.writelines("gauss")

def create_dft_input(mol,tmp,loc):
    # Input Generation by AmberTools
    subprocess.call("conda run -n gmxMMPBSA antechamber -i {0}/{1}_antechamber.mol2 -fi mol2 -o {2}/{1}_gaus.com -gv 1 -ge {1}_gaus.gesp -fo gcrt -pf y -s 0".format(location,mol,tmp),shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    # Change Level of Theory
    correct='''--Link1--
    %chk={0}.chk
    %nproc=14
    %mem=8000MB
    #P wb97xd/6-311G** int(grid=ultrafine) SCF=tight Pop=MK iop(6/33=2) iop(6/42=6) iop(6/50=1) Symmetry=None\n\n'''.format(mol)

    incorrect='''--Link1--\n
            %chk=molecule\n
            #HF/6-31G* SCF=tight Test Pop=MK iop(6/33=2) iop(6/42=6) opt\n
            # iop(6/50=1)\n'''

    new_com = []
    new_com.append(correct)
    with open('{1}/{0}_gaus.com'.format(mol,tmp),'r') as f:
        for i in f.readlines():
            if i in incorrect:
                next
            elif "remark line goes here" in i:
                new_com.append("{0}_gaus.gesp\n\n".format(mol))
            elif "{0}_gaus.gesp".format(mol) in i:
                new_com.append("\n")
                new_com.append(i)
            else:
                new_com.append(i)
    new_com.append(" \n")

    # Write New File
    with open('{1}/{0}_gaus.com'.format(mol,tmp),'w') as f:
        f.writelines(new_com)

def create_complex_topo(mol,tmp,curr):
    # Open Protein and Ligand Structure file

    input_prot = protein
    input_lig = mol

    prot = open('{0}/{1}.gro'.format(curr,input_prot),'r')
    lig = open('{0}/{1}.gro'.format(tmp,input_lig),'r')
    comp = open('{0}/complex.gro'.format(tmp),'w')
    lig_top = open('{0}/{1}.top'.format(tmp,input_lig),'r')
    lig_new_top = open('{0}/{1}.itp'.format(tmp,input_lig),'w')

    s_p = prot.readlines()
    s_l = lig.readlines()

    #Number of atoms
    a_p = int(s_p[1])
    a_l = int(s_l[1])
    a_new = a_p + a_l

    #Coordinates
    c_p = s_p[-1]

    #Header
    header = "Protein Complex for ABFE\n"

    #Structure data
    s_new_p = s_p[2:-1]
    s_new_l = s_l[2:-1]

    #Write New structure file
    comp.write(header)
    comp.write(str(a_new) + "\n")
    comp.writelines(s_new_p)
    comp.writelines(s_new_l)
    comp.write(c_p)

    comp.close()

    #Read and Write LIG topology

    def1 = "[ defaults ]"
    def2 = "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ"
    def3 = "1               2               yes             0.5          0.83333333"

    itp = []

    for i in lig_top.readlines():
        if def1 in i:
            next
        elif def2 in i:
            next
        elif def3 in i:
            next
        else:
            itp.append(i)

    lig_name = itp[-1]

    new_itp = itp[0:-7]

    lig_new_top.writelines(new_itp)
    lig_new_top.close()

    #Read and Write Protein topology


    after ='#include "./amber99sb-ildn_zn.ff/forcefield.itp"'
    include = '\n; Include Ligand Topology\n#include "{}.itp"\n'.format(input_lig)

    n_top = []
    with open('{0}/topol.top'.format(curr),'r') as f:
        for i in f.readlines():
            if after in i:
                n_top.append(i + include)
            else:
                n_top.append(i)
    n_top.append(lig_name)

    with open('{0}/topol.top'.format(tmp),'w') as f:
        f.writelines(n_top)


def dft(mol,tmp,loc,curr):
    # Create Inputs
    create_dft_input(mol,tmp,loc)
    #Run DFT
    subprocess.call("bash {0}/gauss.sh {1}_gaus.com".format(curr,mol),shell=True)
    # RESP charges
    subprocess.call("conda run -n gmxMMPBSA antechamber -i {1}/{0}_gaus.gesp -fi gesp -o {1}/{0}.mol2 -fo mol2 -pf y -s 0 -c resp -eq 2".format(mol,tmp),shell=True, stdout=FNULL, stderr=subprocess.STDOUT)


def convert_pdbqt(location,mol,tmp,form):
    FNULL = open(os.devnull, 'w')
    subprocess.call("obabel {0}/{1}.{2} -O {3}/{1}.pdb -m".format(location,mol,form,tmp),shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    subprocess.call("obabel {1}/{0}1.pdb -O {1}/{0}1.sdf -p 7.4".format(mol,tmp),shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
def conver_for_charge(location,mol,tmp,form):
    FNULL = open(os.devnull, 'w')
    subprocess.call("obabel {0}/{1}.{2} -O {3}/{1}.sdf".format(location,mol,form,tmp),shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

def create_topo(mol,tmp):
    # Parameterize Ligands
    FNULL = open(os.devnull, 'w')
    ligprep_file='''source oldff/leaprc.ff99SB
    source leaprc.gaff 
    loadamberparams {1}/{0}.frcmod 
    lig = loadmol2 {1}/{0}.mol2
    check lig
    saveoff lig {1}/lig.lib
    saveamberparm lig {1}/{0}.prmtop {1}/{0}.rst7
    quit'''.format(mol,tmp) 
    with open("{1}/ligprep{0}.in".format(mol,tmp),"w") as f:
        f.writelines(ligprep_file)
    subprocess.call("parmchk2 -i {1}/{0}.mol2 -f mol2 -o {1}/{0}.frcmod".format(mol,tmp),shell=True,stdout=FNULL,stderr=subprocess.STDOUT)
    #subprocess.call("parmchk2 -i {1}/{0}.mol2 -f mol2 -o {1}/{0}.frcmod".format(mol,tmp),shell=True,stdout=FNULL,stderr=subprocess.STDOUT)
    
    subprocess.call("tleap -s -f {1}/ligprep{0}.in > {1}/ligprep.out".format(mol,tmp),shell=True,stdout=FNULL,stderr=subprocess.STDOUT)

def convert_gmx(mol,tmp):
    # Convert molecule to Gromacs readible file
    par = pmd.load_file("{1}/{0}.prmtop".format(mol,tmp),"{1}/{0}.rst7".format(mol,tmp))
    par.save('{1}/{0}.gro'.format(mol,tmp))
    par.save('{1}/{0}.top'.format(mol,tmp)) 

def bcc(mol,tmp,curr):
    charge = 0
    with open("{1}/{0}1.sdf".format(mol,tmp),"r") as f:
        for i in f.readlines():
            if "M  CHG" in i:
                k = i.split()
                charge = k[-1]
            else:
                next
    os.chdir(tmp)
    m = subprocess.call("antechamber -i {1}/{0}1.sdf -fi sdf -o {1}/{0}.mol2 -fo mol2 -pf y -s 0 -j 5 -c bcc".format(mol,tmp),shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    if m == 1:
        result = subprocess.run("antechamber -i {1}/{0}1.sdf -fi sdf -o {1}/{0}.mol2 -fo mol2 -pf y -s 0 -c bcc -nc {2}".format(mol,tmp,charge),shell=True, capture_output=True)
        try:
            call = result.stdout.decode()
            call2 = result.stderr.decode()
            try:
                os.mkdir(errors)
            except:
                next
            with open("{0}/{1}_err.out".format(errors,mol),"w") as f:
                f.write("###### {0} FAILED #####\nHere is why: \n".format(mol))
                f.write(call)
                f.write(call2)
        except:
            next
    else:
        next
    os.chdir(curr)


def main(i):
    # File formats we are looking for
    mol2 = "mol2"
    pdbqt = "pdbqt"
    moi = i
    conv = moi.replace(".pdbqt","") # JW remove th suffix in the path
    conv = conv.replace(".mol2","")
    conv = conv.replace(".sdf","")


    # Setup TMPDIR
    tmp = os.path.join(os.getcwd(),"tmpdir_{0}".format(conv)) # Make a directory with the mol name
    loc = os.path.join(os.getcwd(),target) # target is path to ligand file
    curr = os.getcwd()
    errors = os.path.join(curr,"errors")
    #checked = os.path.join(curr,"good")
    try:
        shutil.rmtree(tmp)
    except:
        next
    os.mkdir(tmp)
    #os.chdir(tmp)
    
    # Convert PDBQT
    if ".pdbqt" in moi:
        convert_pdbqt(loc,conv,tmp,pdbqt)
    elif ".mol2" in moi:
        if "none" in method:
            shutil.copy("{0}/{1}".format(loc,moi),"{0}/{1}".format(tmp,moi))
        elif "no" in ob:
            shutil.copy("{0}/{1}".format(loc,moi),"{0}/{1}.mol2".format(tmp,conv))
        else:
            convert_pdbqt(loc,conv,tmp,mol2)
    elif ".sdf" in moi:
        shutil.copy("{0}/{1}.sdf".format(loc,conv),"{0}/{1}1.sdf".format(tmp,conv))
    else:
        raise ("ERROR")
    
    moi = conv
    # Parameterize if needed
    if "dft" in method:
        dft(moi,tmp)
    elif "bcc" in method:
        bcc(moi,tmp,curr)
    else:
        next
    
    try:

        # Create Topology
        create_topo(moi,tmp)

        # Convert to Gromacs
        convert_gmx(moi,tmp)
        
        # Create new Destination
        ready = os.path.join(curr,"ready")
        new_d = os.path.join(ready,moi)
        try:
            os.mkdir(new_d)
        except:
            next

        # Create Complex File
        if protein is not None:
            create_complex_topo(moi,tmp,curr)
            shutil.copy("{0}/complex.gro".format(tmp),new_d)
            shutil.copy("{0}/topol.top".format(tmp),new_d)
            shutil.copy("{1}/{0}.itp".format(moi,tmp ),new_d)
        else:
            shutil.copy("complex.gro",new_d)
            shutil.copy("topol.top",new_d)
            next
    except:
        next
    # Finish-up


    os.chdir(curr)
    shutil.rmtree(tmp)


if __name__ == '__main__':
    # Set-up arguments
    global target
    global method
    global protein
    global ob
    global errors

    errors = os.path.join(os.getcwd(),"errors")

    parser = argparse.ArgumentParser(prog='Ligand Parameterizer 1.0',description='This program takes one or multiple ligands and parameterizes them with GAFF using DFT or AM1-BCC charge methods.\n Additionally, it can create complex file if supplied.')
    parser.add_argument('-s','--source')
    parser.add_argument('-m','--method')
    parser.add_argument('-id','--mol_id')
    parser.add_argument('-p','--protein')
    parser.add_argument('-o','--obabel')
    parse = parser.parse_args()

    with_gauss = parse.method
    source = parse.source # JW I assume this is the ligand file
    mol_id = parse.mol_id
    protein = parse.protein
    ob = parse.obabel
   
    target = os.path.join(os.getcwd(),source) # JW target is now path to ligand
    curr = os.getcwd()

    # Determine Method
    if "dft" in with_gauss:
        if "fourier" in socket.gethostname():
            write_gauss_shell(curr)
            method="dft"
            n=4
        else:
            next
    elif "bcc" in with_gauss:
        method="bcc"
        n=None
    elif "none" in with_gauss:
        method="none"
        n=None
    else:
       raise ("Method has to be defined")
    
    # Paralellize process
    print("Starting Checking...")

    if mol_id is None:
        with Pool(n) as pool:
            for i in tqdm(pool.imap(main, os.listdir(target)),total=len(os.listdir(target))):
                next
        #print(str(len(os.listdir("good"))) + "/" + str(len(os.listdir(target))) + " passed.")
        print("Finished Checking.")
    elif mol_id is not None:
        main(mol_id)
        print("Finished Checking.")
    
    