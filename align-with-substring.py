import sys

pos_lines = []
def start(needle):
    for line in sys.stdin:
        if needle not in line: continue

        pos = line.index(needle)

        pos_lines.append((pos, line.strip()))

    if not pos_lines: return
    
    pos_lines.sort(reverse=True)
    max_pos = pos_lines[0][0]
    
    for pos, line in pos_lines:
        print("%s%s" % (
            " "*(max_pos-pos), line))    
        

if __name__ == "__main__":
    start(*sys.argv[1:])
