import sys
import os

from pyLCIO import EVENT, IMPL, IOIMPL, UTIL

try:
    arg1 = sys.argv[3]
except IndexError:
    print "Usage: " + os.path.basename(__file__) + " <inputFile> <outputFile> <maxEvt>"
    sys.exit(1)

# create a reader and open input file
infile = sys.argv[1]
reader=IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(infile)

# create a writer and open the output file
outFile = sys.argv[2]
writer = IOIMPL.LCFactory.getInstance().createLCWriter()
writer.open( outFile, EVENT.LCIO.WRITE_NEW )

maxEvt = sys.argv[3]

i=0
for event in reader:
    i+=1

    outEvent = event

    #print "Max event: " + str(maxEvt)

    if i > int(maxEvt):
        break

    if i%5 == 0:
        print "Saving event " + str(i) + "..."
    writer.writeEvent( outEvent )

reader.close()
writer.flush()
writer.close()
