coverage = config["coverage"]
chromosome = config['chromosome']
seed = config['seed']
folder  = config['folder']

rule all:
    input:
        sdf = directory(f'Data/{folder}/{chromosome}/{chromosome}.sdf'),
        bcf_acc = directory(f'Data/{folder}/benchmarking/accuracy')



rule download_chromosome20:
    output:
        f'Data/{folder}/{chromosome}/{chromosome}.fa.gz'

    shell:
        '''
        mkdir -p Data/{folder}/{chromosome}
        wget \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chromosome}.fa.gz' \
        -O {output}
        '''


rule decompress_reference:
    input:
        f'Data/{folder}/{chromosome}/{chromosome}.fa.gz'

    output:
        f'Data/{folder}/{chromosome}/{chromosome}.fa'

    shell:
        '''
        gunzip -c {input} > {output}
        '''


rule index_reference:
    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa'

    output:
        f'Data/{folder}/{chromosome}/{chromosome}.fa.fai'

    shell:
        '''
        samtools faidx {input.ref}
        '''

rule generate_truth_vcf:
    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa',
        index = f'Data/{folder}/{chromosome}/{chromosome}.fa.fai'

    output:
        vcf = f'Data/{folder}/simulated/truth.vcf'

    shell:
        '''
        mkdir -p Data/{folder}/simulated

        holodeck mutate \
        -r {input.ref} \
        -o {output.vcf} \
        --snp-rate 0.001 \
        --seed {seed}
        '''


rule generate_reads:
    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa',
        vcf = f'Data/{folder}/simulated/truth.vcf'

    output:
        golden_bam = f'Data/{folder}/simulated/sim.golden.bam',
        r1 = f'Data/{folder}/simulated/sim.r1.fastq.gz',
        r2 = f'Data/{folder}/simulated/sim.r2.fastq.gz'

    shell:
        '''
        mkdir -p Data/{folder}/simulated

        holodeck simulate \
        -r {input.ref} \
        -v {input.vcf} \
        -o Data/{folder}/simulated/sim \
        --coverage {coverage} \
        --golden-bam \
        --seed {seed}
        '''


rule build_bowtie_index:
    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa'

    output:
        f'Data/{folder}/bowtie_index/{chromosome}.1.bt2',
        f'Data/{folder}/bowtie_index/{chromosome}.2.bt2',
        f'Data/{folder}/bowtie_index/{chromosome}.3.bt2',
        f'Data/{folder}/bowtie_index/{chromosome}.4.bt2',
        f'Data/{folder}/bowtie_index/{chromosome}.rev.1.bt2',
        f'Data/{folder}/bowtie_index/{chromosome}.rev.2.bt2'

    shell:
        '''
        mkdir -p Data/{folder}/bowtie_index

        bowtie2-build \
        {input.ref} \
        Data/{folder}/bowtie_index/{chromosome}
        '''


rule align_reads:
    input:
        index = f'Data/{folder}/bowtie_index/{chromosome}.1.bt2',
        r1 = f'Data/{folder}/simulated/sim.r1.fastq.gz',
        r2 = f'Data/{folder}/simulated/sim.r2.fastq.gz'

    output:
        bam = f'Data/{folder}/aligned/bowtie2_sorted.bam',
        bam_index = f'Data/{folder}/aligned/bowtie2_sorted.bam.bai'

    shell:
        '''
        mkdir -p Data/{folder}/aligned

        bowtie2 \
        -x Data/{folder}/bowtie_index/{chromosome} \
        -1 {input.r1} \
        -2 {input.r2} \
        -p 4 \
        | samtools sort \
        -o {output.bam} 
        
        samtools index {output.bam}
        '''


rule benchmark_bcftools:
    threads: 1

    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa',
        bam = f'Data/{folder}/aligned/bowtie2_sorted.bam',
        bam_index = f'Data/{folder}/aligned/bowtie2_sorted.bam.bai'
    
    output:
        vcf = f'Data/{folder}/vcfs/bcftools.vcf.gz',
        bench = f'Data/{folder}/benchmarking/bcftools.txt'

    shell:
        '''
        mkdir -p Data/{folder}/vcfs
        mkdir -p Data/{folder}/benchmarking
        /usr/bin/time -v -o {output.bench} \
        bash -c '
            bcftools mpileup \
                -Ou \
                -f {input.ref} \
                {input.bam} \
            | bcftools call \
                -mv \
                -Oz \
                -o {output.vcf}
        '
        '''

rule normalize_vcf:
    input:
        ref = f'Data/{folder}/{chromosome}/{chromosome}.fa',
        truth = f'Data/{folder}/simulated/truth.vcf',
        bcf =  f'Data/{folder}/vcfs/bcftools.vcf.gz'
    
    output:
        norm_truth = f'Data/{folder}/simulated/truth_norm.vcf.gz',
        norm_truth_index = f'Data/{folder}/simulated/truth_norm.vcf.gz.tbi',
        norm_bcf = f'Data/{folder}/vcfs/bcftools_norm.vcf.gz',
        norm_bcf_index = f'Data/{folder}/vcfs/bcftools_norm.vcf.gz.tbi'

    shell:
        '''
        bcftools norm -f {input.ref} -Oz -o {output.norm_truth} {input.truth}
        bcftools index -t {output.norm_truth}
        bcftools norm -f {input.ref} -Oz -o {output.norm_bcf} {input.bcf}
        bcftools index -t {output.norm_bcf}
        '''

rule evaluate_vcf:
    input:
        norm_truth = f'Data/{folder}/simulated/truth_norm.vcf.gz',
        norm_bcf = f'Data/{folder}/vcfs/bcftools_norm.vcf.gz',
        ref =  f'Data/{folder}/{chromosome}/{chromosome}.fa'
    
    output:
        sdf = directory(f'Data/{folder}/{chromosome}/{chromosome}.sdf'),
        output_dir = directory(f'Data/{folder}/benchmarking/accuracy')

    shell:
        '''
        rtg format -o {output.sdf} {input.ref}
        rtg vcfeval \
        -b {input.norm_truth} \
        -c {input.norm_bcf} \
        -t {output.sdf} \
        -o {output.output_dir}

        find {output.output_dir} -type f ! -name 'summary.txt' -delete

        mv {output.output_dir}/summary.txt \
        {output.output_dir}/bcftools_accuracy.txt
        '''
