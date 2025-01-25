import os
import numpy as np

pdb_info = {}
gro_info = {}

for filename in os.listdir('tests/data/fraglib-1r69-fragment_memory_term'):
    if filename[-4:] == '.pdb':
        pdb_info.update({filename[:-4]:[]})
        with open(f'tests/data/fraglib-1r69-fragment_memory_term/{filename}','r') as f:
            coord_lines = []
            for line in f:
                if line[:4] == "ATOM" or line[:6] == "HETATM":
                    x = str(round(float(line[30:38].strip())/10,3))
                    y = str(round(float(line[38:46].strip())/10,3))
                    z = str(round(float(line[46:54].strip())/10,3))
                    pdb_info[filename[:-4]] += [float(x),float(y),float(z)]
                    coord_lines.append(f'{x.strip().ljust(10)}{y.strip().ljust(10)}{z.strip().ljust(10)}\n')
        if not os.path.exists('tests/data/coords_fraglib-1r69-fragment_memory_term'):
            os.mkdir('tests/data/coords_fraglib-1r69-fragment_memory_term')
        #with open(f'tests/data/fraglib-1r69-fragment_memory_term/{filename}.coords','w') as f:
        with open(f'tests/data/coords_fraglib-1r69-fragment_memory_term/{filename}.coords','w') as f:
            for line in coord_lines:
                f.write(line)
    elif filename[-4:] == '.gro':
        gro_info.update({filename[:-4]:[]})
        with open(f'tests/data/fraglib-1r69-fragment_memory_term/{filename}','r') as f:
            coord_lines = []
            for counter,line in enumerate(f):
                if counter >= 2:
                    try:
                        x, y, z, = [coord for coord in line[20:].split(' ') if coord!='']
                    except ValueError:
                        print(filename)
                        print(line)
                        raise
                    gro_info[filename[:-4]] += [float(x.strip()),float(y.strip()),float(z.strip())]
                    coord_lines.append(f'{x.strip().ljust(10)}{y.strip().ljust(10)}{z.strip().ljust(10)}\n')
        if not os.path.exists('tests/data/coords_fraglib-1r69-fragment_memory_term'):
            os.mkdir('tests/data/coords_fraglib-1r69-fragment_memory_term')
        with open(f'tests/data/coords_fraglib-1r69-fragment_memory_term/{filename}.coords','w') as f:
        #with open(f'tests/data/fraglib-1r69-fragment_memory_term/{filename}.coords','w') as f:
            for line in coord_lines:
                f.write(line)

with open('coords_compare.txt','w') as f:
    combined_info = {}
    failed = []
    for pdbid in pdb_info.keys():
        try:
            combined_info[pdbid] = np.array(pdb_info[pdbid]) - np.array(gro_info[pdbid])
        except ValueError:
            failed.append(pdbid)
            continue
        max_num = np.max(np.abs(combined_info[pdbid]))
        if max_num > 0.0010000001:
            f.write(f'{pdbid}\n')
            f.write(f'max: {max_num}\n')
        #f.write(f'min: {np.min(np.abs(combined_info[pdbid]))}\n')
        if np.min(np.abs(combined_info[pdbid])) != 0:
            raise ValueError(f"min was {np.min(np.abs(combined_info[pdbid]))} for pdbid {pdbid}!")
    f.write('failed\n')
    f.write(str(failed))
    if failed:
        raise AssertionError("strange issue in coordinate comparison. see coords_compare.txt")
