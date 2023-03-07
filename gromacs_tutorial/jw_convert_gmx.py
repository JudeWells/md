import argparse
import parmed as pmd

"""
JW script that factors out the convert gmx python part 
of Daniel's script parameterize_molz.py
"""

def convert_gmx(inpath):
    # Convert molecule to Gromacs readible file
    par = pmd.load_file(inpath, inpath.replace('.prmtop', '.rst7'))
    par.save(inpath.replace('.prmtop', '.gro'))
    par.save(inpath.replace('.prmtop', '.top'))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='path to .prmtop file')
    args = parser.parse_args()
    convert_gmx(args.input)

