import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored
pymol.finish_launching()

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 4)   # Lowest resolution

pdb = '1y57_human_numbering'

cmd.load('../pdbs/' + pdb + '.pdb')  # use the name of your pdb file
stored.residues = []
cmd.iterate('name ca', 'stored.residues.append(resi)')

sasa_per_residue = []
idx = []
for i in stored.residues:
    sasa_per_residue.append(cmd.get_area('resi %s' % i))
    idx.append(i)

with open(pdb + '_SASA.txt', 'w') as f:
    for i, item in zip(idx, sasa_per_residue):
        #f.write("%s\n" % item)
        line = "{}\t{}\n".format(i, item)
        f.write(line)

print sum(sasa_per_residue)
print cmd.get_area('all')  # just to check that the sum of sasa per residue equals the total area
