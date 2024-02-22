import pyfastx
from threading import Thread
fq = pyfastx.Fastq('../test/ecoli.fq.gz')
qc_all_read = dict()
def get_read_number( fq, qc_all_read ):
    qc_all_read['read_number'] = len(fq)
    return qc_all_read
def get_base_number( fq, qc_all_read ):
    qc_all_read['base_number'] = fq.size
    return qc_all_read
def get_gc_content( fq, qc_all_read ):
    qc_all_read['gc_content'] = fq.gc_content
    return qc_all_read
def get_base_composition( fq, qc_all_read ):
    qc_all_read['base_composition'] = fq.composition
    return qc_all_read
def get_read_minimum_length( fq, qc_all_read ):
    qc_all_read['read_minimum_length'] = fq.minlen
    return qc_all_read
def get_read_maximum_length( fq, qc_all_read ):
    qc_all_read['read_maximum_length'] = fq.maxlen
    return qc_all_read
def get_base_minimum_quality( fq, qc_all_read ):
    qc_all_read['base_minimum_quality'] = fq.minqual - fq.phred
    return qc_all_read
def get_base_maximum_quality( fq, qc_all_read ):
    qc_all_read['base_maximum_quality'] = fq.maxqual - fq.phred
    return qc_all_read
def quick_qc_table(fq, qc_all_read):
    threads = []
    threads.append(Thread(target=get_read_number, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_base_number, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_gc_content, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_base_composition, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_read_minimum_length, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_read_maximum_length, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_base_minimum_quality, args=(fq, qc_all_read)))
    threads.append(Thread(target=get_base_maximum_quality, args=(fq, qc_all_read)))
    for fun in threads:
        fun.start()
    return qc_all_read

qc_all_read_table = quick_qc_table(fq, qc_all_read)

for k, v in qc_all_read_table.items():
    print(k, v)
