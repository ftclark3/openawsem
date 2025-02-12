import subprocess
import os

# read in names of fragments that should be present in our test
fraglist_file = '/home/fc36/openawsem/tests/data/1r69-fragment_memory_term.mem'
names = []
with open(fraglist_file,'r') as f:
    for counter, line in enumerate(f):
        if counter < 4:
            continue
        else:
            names.append(line.split(' ')[0].split('/')[-1].split('.')[0])

# ensure that our fraglib directory contains a pdb for each of the names (and only the names) in the list
for other_filename in os.listdir('/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/'):
    other_name = other_filename.split('.')[0]
    if other_name in names:
        subprocess.run(['cp',f'cif_handling/carlos-master_create/process_cif/newest/1r69/fraglib/{other_name}.pdb',
            f'/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/{other_name}.pdb'])

new_names = [filename.split('.')[0] for filename in os.listdir('/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/')]
for new_name in new_names:
    pdb = os.path.isfile(f'/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/{new_name}.pdb')
    gro = os.path.isfile(f'/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/{new_name}.gro')
    if pdb and gro:
        pass
    elif pdb and not gro:
        subprocess.run(['rm',f'/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/{new_name}.pdb'])
    elif gro and not pdb:
        subprocess.run(['rm',f'/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/{new_name}.gro'])
    else:
        raise AssertionError(new_name) 

new_names = [filename.split('.')[0] for filename in os.listdir('/home/fc36/openawsem/tests/data/fraglib-1r69-fragment_memory_term/')]
lines = []
hits = []
with open('/home/fc36/openawsem/tests/data/1r69-fragment_memory_term.mem','r') as f:
    for counter, line in enumerate(f):
        if counter < 4:
            lines.append(line)
        else:
            found_flag = False
            for name in new_names:
                if name in line:
                    found_flag = True
                    if name not in hits:
                        lines.append(line)
                    hits.append(name)
     
with open('/home/fc36/openawsem/tests/data/1r69-fragment_memory_term.mem_foo','w') as f:
    for line in lines:
        f.write(line)