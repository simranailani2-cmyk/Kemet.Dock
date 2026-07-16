def parse_vina_output(output_str):
    lines = output_str.split('\n')
    data = []
    parsing = False
    for line in lines:
        if line.startswith('-----+------------+----------+----------'):
            parsing = True
            continue
        if parsing:
            parts = line.split()
            if len(parts) == 4 and parts[0].isdigit():
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    # user asked for 'RMSD', we can use rmsd u.b.
                    rmsd_ub = float(parts[3])
                    data.append({
                        'Binding Mode': mode,
                        'Affinity (kcal/mol)': affinity,
                        'RMSD': rmsd_ub
                    })
                except ValueError:
                    pass
            elif len(parts) == 0 or 'Writing' in line:
                break
    return data

sample_out = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1         -7.1      0.000      0.000
   2         -6.8      2.322      5.144
"""
print(parse_vina_output(sample_out))
