import os, sys

# if __name__ == '__main__':
print(sys.argv)
with open(name = sys.argv[1], mode="r") as fin:
    with open(name = sys.argv[2], mode="w+") as fout:
        for line in fin:
	    print 'line',line;	
	    filledBucket = line.split(' ')[2];
            print 'filled',filledBucket
            valid = ' 1.0 1.0 1.0 1.0\n';
            invalid = ' 0.0 0.0 0.0 0.0\n'
            lineOut = line.rstrip() + valid;
            if filledBucket=="0":
                lineOut=line.rstrip() + invalid; 	
	    print lineOut;	
            fout.write(lineOut);
