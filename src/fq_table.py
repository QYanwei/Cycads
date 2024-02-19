import pyfastx

fq = pyfastx.Fastq('../test/ecoli.fq.gz')

print('read number:', len(fq))
print('base number:', fq.size)
print('gc content:', fq.gc_content)
print('base composition:', fq.composition)
print('average length:', fq.avglen)
print('maximum length:', fq.maxlen)
print('minimum length:', fq.minlen)
print('base minimum quality:', fq.minqual - fq.phred)
print('base maximum quality:', fq.maxqual - fq.phred)


