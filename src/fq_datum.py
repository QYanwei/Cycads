import pyfastx

fq = pyfastx.Fastx('../test/ecoli.fq.gz')

i = 0
for name, seq, qual in fq:
    print( name )
    print( len(seq) )
    print( len(qual) )
    i += 1

print(i)