rule all:
    input:
        'Data/aligned/bowtie2_sorted.bam',
        'Data/aligned/bowtie2_sorted.bam.bai'


rule download_chromosome20:
    output:
        'Data/chr20/chr20.fa.gz'

    shell:
        '''
        mkdir -p Data/chr20
        wget \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz' \
        -O {output}
        '''


rule decompress_reference:
    input:
        'Data/chr20/chr20.fa.gz'

    output:
        'Data/chr20/chr20.fa'

    shell:
        '''
        gunzip -c {input} > {output}
        '''


rule index_reference:
    input:
        ref = 'Data/chr20/chr20.fa'

    output:
        'Data/chr20/chr20.fa.fai'

    shell:
        '''
        samtools faidx {input.ref}
        '''

rule generate_truth_vcf:
    input:
        ref = 'Data/chr20/chr20.fa',
        index = 'Data/chr20/chr20.fa.fai'

    output:
        vcf = 'Data/simulated/truth.vcf'

    shell:
        '''
        mkdir -p Data/simulated

        holodeck mutate \
        -r {input.ref} \
        -o {output.vcf} \
        --snp-rate 0.001
        '''


rule generate_reads:
    input:
        ref = 'Data/chr20/chr20.fa',
        vcf = 'Data/simulated/truth.vcf'

    output:
        golden_bam = 'Data/simulated/sim.golden.bam',
        r1 = 'Data/simulated/sim.r1.fastq.gz',
        r2 = 'Data/simulated/sim.r2.fastq.gz'

    shell:
        '''
        mkdir -p Data/simulated

        holodeck simulate \
        -r {input.ref} \
        -v {input.vcf} \
        -o Data/simulated/sim \
        --coverage 10 \
        --golden-bam
        '''


rule build_bowtie_index:
    input:
        ref = 'Data/chr20/chr20.fa'

    output:
        'Data/bowtie_index/chr20.1.bt2',
        'Data/bowtie_index/chr20.2.bt2',
        'Data/bowtie_index/chr20.3.bt2',
        'Data/bowtie_index/chr20.4.bt2',
        'Data/bowtie_index/chr20.rev.1.bt2',
        'Data/bowtie_index/chr20.rev.2.bt2'

    shell:
        '''
        mkdir -p Data/bowtie_index

        bowtie2-build \
        {input.ref} \
        Data/bowtie_index/chr20
        '''


rule align_reads:
    input:
        index = 'Data/bowtie_index/chr20.1.bt2',
        r1 = 'Data/simulated/sim.r1.fastq.gz',
        r2 = 'Data/simulated/sim.r2.fastq.gz'

    output:
        bam = 'Data/aligned/bowtie2_sorted.bam',
        bam_index = 'Data/aligned/bowtie2_sorted.bam.bai'

    shell:
        '''
        mkdir -p Data/aligned

        bowtie2 \
        -x Data/bowtie_index/chr20 \
        -1 {input.r1} \
        -2 {input.r2} \
        -p 4 \
        | samtools sort \
        -o {output.bam} 
        
        samtools index {output.bam}
        '''



