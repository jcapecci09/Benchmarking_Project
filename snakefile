coverage = config["coverage"]
chromosome = config['chromosome']
seed = config['seed']
folder  = config['folder']

rule all:
    input:
        f'Data/{chromosome}/{folder}/benchmarking/collected.done'


rule download_chromosome20:
    output:
        f'Data/{chromosome}/{chromosome}.fa.gz'

    shell:
        '''
        mkdir -p Data/{chromosome}
        wget \
        'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chromosome}.fa.gz' \
        -O {output}
        '''


rule decompress_reference:
    input:
        f'Data/{chromosome}/{chromosome}.fa.gz'

    output:
        f'Data/{chromosome}/{chromosome}.fa'

    shell:
        '''
        gunzip -c {input} > {output}
        '''


rule index_reference:
    input:
        ref = f'Data/{chromosome}/{chromosome}.fa'

    output:
        f'Data/{chromosome}/{chromosome}.fa.fai'

    shell:
        '''
        samtools faidx {input.ref}
        '''

rule generate_truth_vcf:
    input:
        ref = f'Data/{chromosome}/{chromosome}.fa',
        index = f'Data/{chromosome}/{chromosome}.fa.fai'

    output:
        vcf = f'Data/{chromosome}/{folder}/simulated/truth.vcf'

    shell:
        '''
        mkdir -p Data/{chromosome}/{folder}/simulated

        holodeck mutate \
        -r {input.ref} \
        -o {output.vcf} \
        --snp-rate 0.001 \
        --seed {seed}
        '''


rule generate_reads:
    input:
        ref = f'Data/{chromosome}/{chromosome}.fa',
        vcf = f'Data/{chromosome}/{folder}/simulated/truth.vcf'

    output:
        golden_bam = f'Data/{chromosome}/{folder}/simulated/sim.golden.bam',
        r1 = f'Data/{chromosome}/{folder}/simulated/sim.r1.fastq.gz',
        r2 = f'Data/{chromosome}/{folder}/simulated/sim.r2.fastq.gz'

    shell:
        '''
        mkdir -p Data/{chromosome}/{folder}/simulated

        holodeck simulate \
        -r {input.ref} \
        -v {input.vcf} \
        -o Data/{chromosome}/{folder}/simulated/sim \
        --coverage {coverage} \
        --golden-bam \
        --seed {seed}
        '''


rule build_bowtie_index:
    input:
        ref = f'Data/{chromosome}/{chromosome}.fa'

    output:
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.1.bt2',
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.2.bt2',
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.3.bt2',
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.4.bt2',
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.rev.1.bt2',
        f'Data/{chromosome}/{folder}/bowtie_index/{chromosome}.rev.2.bt2'

    shell:
        '''
        mkdir -p Data/{chromosome}/{folder}/bowtie_index

        bowtie2-build \
        {input.ref} \
        Data/{chromosome}/{folder}/bowtie_index/{chromosome}
        '''


rule align_reads:
    input:
        ref = f'Data/{chromosome}/{chromosome}.fa',
        r1 = f'Data/{chromosome}/{folder}/simulated/sim.r1.fastq.gz',
        r2 = f'Data/{chromosome}/{folder}/simulated/sim.r2.fastq.gz'

    output:
        bam = f'Data/{chromosome}/{folder}/aligned/bowtie2_sorted.bam',
        bam_index = f'Data/{chromosome}/{folder}/aligned/bowtie2_sorted.bam.bai'

    shell:
        '''
        mkdir -p Data/{chromosome}/{folder}/aligned

        minimap2 \
            -ax sr \
            {input.ref} \\
            {input.r1} \
            {input.r2} \
            -t 4 \
            | samtools sort \
            -o {output.bam}

        samtools index {output.bam}
        '''


rule benchmark_bcftools:
    threads: 1

    input:
        ref = f'Data/{chromosome}/{chromosome}.fa',
        bam = f'Data/{chromosome}/{folder}/aligned/bowtie2_sorted.bam',
        bam_index = f'Data/{chromosome}/{folder}/aligned/bowtie2_sorted.bam.bai'
    
    output:
        vcf = f'Data/{chromosome}/{folder}/vcfs/bcftools.vcf.gz',
        bench = f'Data/{chromosome}/{folder}/benchmarking/bcftools.txt'

    shell:
        '''
        mkdir -p Data/{chromosome}/{folder}/vcfs
        mkdir -p Data/{chromosome}/{folder}/benchmarking
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
        ref = f'Data/{chromosome}/{chromosome}.fa',
        truth = f'Data/{chromosome}/{folder}/simulated/truth.vcf',
        bcf =  f'Data/{chromosome}/{folder}/vcfs/bcftools.vcf.gz'
    
    output:
        norm_truth = f'Data/{chromosome}/{folder}/simulated/truth_norm.vcf.gz',
        norm_truth_index = f'Data/{chromosome}/{folder}/simulated/truth_norm.vcf.gz.tbi',
        norm_bcf = f'Data/{chromosome}/{folder}/vcfs/bcftools_norm.vcf.gz',
        norm_bcf_index = f'Data/{chromosome}/{folder}/vcfs/bcftools_norm.vcf.gz.tbi'

    shell:
        '''
        bcftools norm -f {input.ref} -Oz -o {output.norm_truth} {input.truth}
        bcftools index -t {output.norm_truth}
        bcftools norm -f {input.ref} -Oz -o {output.norm_bcf} {input.bcf}
        bcftools index -t {output.norm_bcf}
        '''


rule evaluate_vcf:
    input:
        norm_truth = f'Data/{chromosome}/{folder}/simulated/truth_norm.vcf.gz',
        norm_bcf = f'Data/{chromosome}/{folder}/vcfs/bcftools_norm.vcf.gz',
        ref =  f'Data/{chromosome}/{chromosome}.fa'
    
    output:
        sdf = directory(f'Data/{chromosome}/{folder}/{chromosome}.sdf'),
        output_dir = directory(f'Data/{chromosome}/{folder}/benchmarking/accuracy')

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


rule append_to_csv:
    input:
        bench = f'Data/{chromosome}/{folder}/benchmarking/bcftools.txt',
        acc_dir = f'Data/{chromosome}/{folder}/benchmarking/accuracy'

    params:
        bench_dir = f'Data/{chromosome}/{folder}/benchmarking'
    
    output:
        f'Data/{chromosome}/{folder}/benchmarking/collected.done'
    
    shell:
        '''
        python scripts/collect_stats.py -i {params.bench_dir} -r {chromosome} -c {coverage} -s {seed}
        touch {output}
        '''
