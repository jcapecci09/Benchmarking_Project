rule all:
    input:
        'Data/simulated/sim.golden.bam',
        'Data/simulated/sim.r1.fastq.gz',
        'Data/simulated/sim.r2.fastq.gz',
        'Data/aligned/sim.sorted.bam'


rule download_chromosome22:
    output:
        'Data/chr22/chr22.fa.gz'

    shell:
        '''
        mkdir -p Data/chr22
        wget --timestamping \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz' \
        -O {output}
        '''


rule decompress_reference:
    input:
        'Data/chr22/chr22.fa.gz'

    output:
        'Data/chr22/chr22.fa'

    shell:
        '''
        gunzip -c {input} > {output}
        '''


rule index_reference:
    input:
        ref = 'Data/chr22/chr22.fa'

    output:
        'Data/chr22/chr22.fa.fai'

    shell:
        '''
        samtools faidx {input.ref}
        '''

rule generate_truth_vcf:
    input:
        ref = 'Data/chr22/chr22.fa',
        index = 'Data/chr22/chr22.fa.fai'

    output:
        vcf = 'Data/simulated/truth.vcf'

    shell:
        '''
        mkdir -p Data/simulated

        ./holodeck/target/release/holodeck mutate \
        -r {input.ref} \
        -o {output.vcf} \
        --snp-rate 0.001
        '''


rule generate_reads:
    input:
        ref = 'Data/chr22/chr22.fa',
        vcf = 'Data/simulated/truth.vcf'

    output:
        golden_bam = 'Data/simulated/sim.golden.bam',
        r1 = 'Data/simulated/sim.r1.fastq.gz',
        r2 = 'Data/simulated/sim.r2.fastq.gz'

    shell:
        '''
        mkdir -p Data/simulated

        ./holodeck/target/release/holodeck simulate \
        -r {input.ref} \
        -v {input.vcf} \
        -o Data/simulated/sim \
        --coverage 10 \
        --golden-bam
        '''


rule build_bowtie_index:
    input:
        ref = 'Data/chr22/chr22.fa'

    output:
        'Data/bowtie_index/chr22.1.bt2',
        'Data/bowtie_index/chr22.2.bt2',
        'Data/bowtie_index/chr22.3.bt2',
        'Data/bowtie_index/chr22.4.bt2',
        'Data/bowtie_index/chr22.rev.1.bt2',
        'Data/bowtie_index/chr22.rev.2.bt2'

    shell:
        '''
        mkdir -p Data/bowtie_index

        bowtie2-build \
        {input.ref} \
        Data/bowtie_index/chr22
        '''


rule align_reads:
    input:
        index = 'Data/bowtie_index/chr22.1.bt2',
        r1 = 'Data/simulated/sim.r1.fastq.gz',
        r2 = 'Data/simulated/sim.r2.fastq.gz'

    output:
        bam = 'Data/aligned/sim.sorted.bam'

    shell:
        '''
        mkdir -p Data/aligned

        bowtie2 \
        -x Data/bowtie_index/chr22 \
        -1 {input.r1} \
        -2 {input.r2} \
        -p 4 \
        | samtools sort \
        -o {output.bam} \
        
        samtools index {output.bam}
        '''