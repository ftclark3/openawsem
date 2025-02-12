lines = []
with open('4tlt/frags.mem','r') as f:
    for line in f:
        if 'fraglib' in line:
            line_split = line.split('fraglib')
            assert len(line_split) == 2, line_split
            line = f'{line_split[0]}fraglib-4tlt-fragment_memory_term{line_split[1]}'
        lines.append(line)

with open('4tlt/4tlt-fragment_memory_term.mem','w') as f:
    for line in lines:
        f.write(line)
 
