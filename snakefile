rule all:
    input:
        'Data/sim.golden.bam',
        'Data/sim.r1.fastq.gz',
        'Data/sim.r2.fastq.gz'


rule download_chromosome22:
    output:
        'Data/chr22.fa.gz'

    shell:
        '''
        mkdir -p Data
        wget --timestamping \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz' \
        -O {output}
        '''


rule decompress_reference:
    input:
        'Data/chr22.fa.gz'

    output:
        'Data/chr22.fa'

    shell:
        '''
        gunzip -c {input} > {output}
        '''


rule index_reference:
    input:
        ref = 'Data/chr22.fa'

    output:
        'Data/chr22.fa.fai'

    shell:
        '''
        samtools faidx {input.ref}
        '''


rule generate_truth_vcf:
    input:
        ref = 'Data/chr22.fa'

    output:
        vcf = 'Data/mutations.vcf'

    shell:
        '''
        ./holodeck/target/release/holodeck mutate \
        -r {input.ref} \
        -o {output.vcf} \
        --snp-rate 0.001
        '''

rule generate_reads:
    input:
        vcf = 'Data/mutations.vcf',
        ref = 'Data/chr22.fa'
    
    output:
        'Data/sim.golden.bam',
        'Data/sim.r1.fastq.gz',
        'Data/sim.r2.fastq.gz'
    
    shell:
        '''
        ./holodeck/target/release/holodeck simulate \
        -r {input.ref} \
        -v {input.vcf} \
        -o Data/sim \
        --coverage 30 \
        --golden-bam
        '''