import os, sys

# if __name__ == '__main__':

with open(name = sys.argv[1], mode="r") as fin:
    with open(name = sys.argv[2], mode="w+") as fout:
        for line in fin:
            if line.startswith('*'):
                line = 'c' + line[1:] 
            fout.write(line);
