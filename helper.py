import re

def find_session_variables():
    with open('app.py', 'r') as f:
        content = f.read()

    # Let's extract the loop over phytochemicals to see what variables are available
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if "for idx, phyto_name in enumerate" in line:
            start_idx = i
            break

    print("Found loop start at line:", start_idx)

    # We want to see how receptor name and plant name are accessed
    for i in range(start_idx, min(start_idx + 100, len(lines))):
        if "phyto_name" in lines[i] or "receptor" in lines[i].lower():
            print(f"{i}: {lines[i]}")

if __name__ == '__main__':
    find_session_variables()
