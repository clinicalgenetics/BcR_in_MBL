import os
import glob
import pysam
import snakemake
import subprocess as sp
import pandas as pd
import numpy as np



print('start')
cwd = os.getcwd()
input_dir = os.path.join("/path/to/cram")
out_path = os.path.join('/path/to/output', 'pipeline_output')
conditions = [os.path.basename(i) for i in glob.glob(os.path.join(input_dir, "*-MBL-*"))]
samples = [j.split("/")[-2] for i in conditions for j in glob.glob(os.path.join(input_dir, i, "*/*.cram"))]
data_dict = dict(zip(conditions, samples))
print(samples)
print(conditions)
print(data_dict)
config = {"data":data_dict,
          "indir": input_dir,
          "outdir": out_path,
      	  "ref": "/path/to/human_g1k_v37.fasta",
          "panel_regions": "/path/to/gmslymphoid_7.2_hg19_design.bed",
          "chrs_size": "/path/to/hg19.chrom.sizes",
          "bwa_idx": "/path/to/bwa_idx/bwa",
          "cll_genes": "/path/to/cllgenes_name_overlap_gms.tsv"}


rule all:
    input:
       [expand(os.path.join(config['outdir'], "bam_rx", "{condition}_{sample}_RX.bam"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "bam_tagged", "{condition}_{sample}_tagged.bam"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.bam"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_R2_consensus.fastq.gz"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.realign.bam"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "vcf_files", "{condition}_{sample}.vcf"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.reheader.vcf"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.vcf"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.filtered.vcf"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.filtered.tsv"),
       condition=key, sample=value) for key, value in data_dict.items()],
       [expand(os.path.join(config['outdir'], "{condition}_{sample}.hs_metrics.txt"),
       condition=key, sample=value) for key, value in data_dict.items()],


rule cram2bam:
    input:
        os.path.join(config['indir'], "{condition}", "{sample}", "tumor.merged.cram")
    output:
        os.path.join(config['outdir'], "{condition}_{sample}.bam")
    params:
        config["ref"]
    shell:
        """
        samtools view -@ 8 -b -T {params} -o {output} {input}; samtools index {output} 
        """


rule add_RX:
    input:
        os.path.join(config['outdir'], "{condition}_{sample}.bam")
    output:
        os.path.join(config['outdir'], "bam_rx", "{condition}_{sample}_RX.bam")
    run:
        bam = pysam.AlignmentFile(input[0], 'rb')
        out_bam = pysam.AlignmentFile(output[0], "wb", template=bam)
        iter = bam.fetch(until_eof=True, multiple_iterators=True)
        for rec in iter:
            fixed_umi = rec.qname.split(":")[-1]
            fixed_umi = fixed_umi.split("_")[1:]
            A = fixed_umi[0][0:3]
            B = fixed_umi[1][0:3]
            fixed_umi = '-'.join([A,B])            
	    rec.set_tag("RX", fixed_umi)
            out_bam.write(rec)
        bam.close()
        out_bam.close()


rule sort_tag:
    input:
        rx_bam = os.path.join(config['outdir'], "bam_rx", "{condition}_{sample}_RX.bam")
    output:
        tagged = os.path.join(config['outdir'], "bam_tagged", "{condition}_{sample}_tagged.bam")
    shell:
        """
        fgbio SortBam -s Queryname -i {input} -o /dev/stdout |\
        fgbio SetMateInformation -i /dev/stdin --allow-missing-mates false -o {output.tagged}
        """

#TODO add thread!
#Diff with old version: --min-map-q 20 instead of 1, --allow-inter-contig false  instead of true
rule umi_group:
    input:
        os.path.join(config['outdir'], "bam_tagged", "{condition}_{sample}_tagged.bam")
    output:
        hist = os.path.join(config['outdir'], "bam_umigroup", "{condition}_{sample}_hist.size.paired"),
        grouped = os.path.join(config['outdir'], "bam_umigroup", "{condition}_{sample}_umigrouped.bam")
    threads:
        8
    shell:
        """
        fgbio GroupReadsByUmi -i {input} -s paired --edits 1 --min-map-q 20 \
        --allow-inter-contig false --include-non-pf-reads false  \
        --family-size-histogram {output.hist} -o {output.grouped}
        """

# Diff --error-rate-post-umi 30  instead of 40, --min-reads 1 0 0 instead of  3 1 1
rule make_consensus:
    input:
        os.path.join(config['outdir'], "bam_umigroup", "{condition}_{sample}_umigrouped.bam")
    output:
        bam = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.bam")
    threads:
        8
    shell:
        """
        fgbio SortBam -s TemplateCoordinate -i {input} -o /dev/stdout | \
        fgbio CallDuplexConsensusReads --min-reads 1 0 0 --threads {threads} --error-rate-pre-umi 45 \
        --error-rate-post-umi 30  --trim false -i /dev/stdin --min-input-base-quality 20 \
        --consensus-call-overlapping-bases true  -o {output.bam}
        """

# Diff -M 1 instead of 4 , -E 0.025 instead 0.010
rule filter_consensus:
    input:
        os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.bam")
    output:
        R1_consensus = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_R1_consensus.fastq.gz"),
        R2_consensus = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_R2_consensus.fastq.gz")
    params:
        ref = config["ref"],
        tmp_dir ="/tmp/"
    threads:
        8
    shell:
        """
        fgbio -Xmx16g FilterConsensusReads -i {input} -o /dev/stdout -r {params.ref} -M 1 -N 1 -E 0.025 \
        | samtools sort -@ {threads} -n \
        | samtools fastq -1 {output.R1_consensus} -2 {output.R2_consensus} -
        """

rule align_consensus:
    input:
        R1_consensus = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_R1_consensus.fastq.gz"),
        R2_consensus = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_R2_consensus.fastq.gz")
    output:
        bam = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.realign.bam"),
        index = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.realign.bam.bai")
    threads:
        8
    params:
        ref = config["bwa_idx"]
    shell:
        """
        bwa-mem2 mem -t {threads} {params.ref} -M {input.R1_consensus} {input.R2_consensus}\
        | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


rule read_panel:
    output:
        temp("tmp.bed")
    params:
        panel = config['panel_regions'],
        chrs_size = config['chrs_size']

    run:
        gms_lym = pd.read_csv(params['panel'], sep = '\t', header = None)
        chrs = np.unique(gms_lym[0].values)
        chr_size = pd.read_csv(params['chrs_size'], sep = '\t', header = None)
        chr_path = os.path.dirname(params['chrs_size'])
        df2save = pd.DataFrame()
        for ch in chrs:
             print(ch)
             length = chr_size.loc[chr_size[0] == "chr"+ch][1].values[0]
             print(length)
             gms_chr = gms_lym.loc[gms_lym[0] == ch]
             gms_chr[1] = np.maximum(gms_chr[1].astype(int)-100, 0)
             gms_chr[2] = np.minimum(gms_chr[2].astype(int)+100, length-1)
             if len(df2save) == 0:
                df2save=gms_chr
             else:
                df2save = pd.concat([df2save, gms_chr])

        df2save.to_csv(output[0], sep = "\t", header = None, index = False)

rule sort_panel:
    input:
        tmp = "tmp.bed"
    output:
        regions = "region.bed"
    shell:
        """
        bedtools sort -i {input.tmp}| bedtools merge > {output.regions}
        """


# Diff I was 50 instead of 200.
rule run_vardict:
    input:
        bam = os.path.join(config['outdir'], "bam_consensus", "{condition}_{sample}_consensus.realign.bam"),
	regions = "region.bed"
    threads:
        8
    params:
        ref = config["ref"],
        tmp_dir ="/tmp/"
    output:
        vcf_file = os.path.join(config['outdir'], "vcf_files", "{condition}_{sample}.vcf"),
    shell:
        """
        VarDict -I 200 -G  {params.ref}  -f 0.01 -N TUMOR \
        -th {threads}  \
        -b  {input.bam} -c 1 -S 2 -E 3 {input.regions}  | \
        teststrandbias.R | var2vcf_valid.pl -N TUMOR -E -f 0.01 > \
        {output.vcf_file}
        """



rule vcf_reheader:
    input:
        vcf = os.path.join(config['outdir'], "vcf_files", "{condition}_{sample}.vcf")
    params:
        fai = config['ref']+".fai"
    output:
        vcf_reheaded = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.reheader.vcf")
    shell:
        """
        bcftools reheader --fai {params.fai} {input.vcf} | bcftools view -f PASS | bcftools sort - -o {output.vcf_reheaded}
        """


rule run_vep:
    input:
        vcf_file = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.reheader.vcf")
    params:
        ref_fasta = config['ref'] 
    output:
        annot_vcf = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.vcf"),
        html = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotation.summary.html")
    shell:
        """
        vep -i {input.vcf_file} --format 'vcf' \
        --use_transcript_ref \
        --hgvsg  \
        --offline \
        --verbose \
        --cache \
        --force_overwrite  \
        --vcf \
        --assembly GRCh37 \
        --allow_non_variant \
        --dont_skip \
        --refseq \
        --output_file {output.annot_vcf} \
        --stats_file {output.html} \
        --species homo_sapiens \
        --dir_cache $VEP_CACHE \
        --everything \
        --individual_zyg all \
        --sift b  --polyphen b \
        --fasta {params.ref_fasta} \
        --plugin CADD,$VEP_CACHE/CADD/GRCh37/whole_genome_SNVs.tsv.gz,$VEP_CACHE/CADD/GRCh37/gnomad.genomes.r2.1.1.indel.tsv.gz \
        --plugin AlphaMissense,file=$VEP_CACHE/AlphaMissense/AlphaMissense_hg19.tsv.gz \
        --synonyms $VEP_CACHE/homo_sapiens/111_GRCh37/chr_synonyms.txt
        """


rule filter_vcf:
      input:
         vcf = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.vcf")
      params:
         cll_genes = config['cll_genes']
      output:
         filtered_vcf = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.filtered.vcf")
      shell:
         """
         filter_vep --only_matched --filter  \"canonical and SYMBOL in {params.cll_genes}\" \
         --input_file {input.vcf} \
         | bcftools view --apply-filter PASS \
         | bcftools filter --exclude \"INFO/AF<0.01\" \
         | bcftools view --types snps,mnps,indels --output-type v --output-file {output.filtered_vcf}
         """



rule split_vep:
      input:
         filtered_vcf = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.filtered.vcf")
      output:
         filtered_tsv = os.path.join(config['outdir'], "annotated_vcfs", "{condition}_{sample}.annotated.filtered.tsv")
      shell:
         """
         bcftools +split-vep {input.filtered_vcf} -d -A tab -H -f \
         '%SAMPLE\t%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%TYPE\t'\
         '%DP\t%VD\t%BIAS\t%REFBIAS\t%VARBIAS\t%CSQ\t[%GT]\t[%DP]\t'\
         '[%VD]\t[%AD]\t[%AF]\t[%RD]\t[%ALD]\n' \
         > {output.filtered_tsv}
         """

rule make_refdict:
      input:
         bed = "region.bed"
      params:
         tmp_dir ="/tmp/",
         mem = '16g',
         fa = config['ref']
      output:
         ref_dict = os.path.join(config['outdir'], "reference.dict"),
      shell:
         """
         module load bioinfo-tools picard/3.1.1;
         java -Xmx{params.mem} -jar $PICARD  \
         CreateSequenceDictionary \
         -R {params.fa} \
         -O {output.ref_dict}
         """

rule picard_bedinterval:
      input:
         bed = "region.bed",
         ref_dict = os.path.join(config['outdir'], "reference.dict"),
      params:
         tmp_dir ="/tmp/",
         mem = '16g',
      output:
         bedintervals = os.path.join(config['outdir'], "picard.bedintervals")
      shell:
         """
         module load bioinfo-tools picard/3.1.1;
         java -Xmx{params.mem} -jar $PICARD \
         BedToIntervalList \
         I={input.bed} \
         O={output.bedintervals} \
         SD={input.ref_dict};
         """


rule run_picard:
      input:
         bam = os.path.join(config['outdir'], "{condition}_{sample}.bam"),
         bedintervals = os.path.join(config['outdir'], "picard.bedintervals")
      params:
         tmp_dir ="/tmp/",
         mem = '16g',
         fa = config['ref']
      output:
         hs_metrics = os.path.join(config['outdir'], "{condition}_{sample}.hs_metrics.txt")
      shell:
         """
         module load bioinfo-tools picard/3.1.1;
         java -Djava.io.tmpdir={params.tmp_dir} -Xmx{params.mem} -jar $PICARD \
         CollectHsMetrics \
         BI={input.bedintervals} \
         TI={input.bedintervals} \
         I={input.bam} \
         O={output.hs_metrics} \
         R={params.fa} \
         COVERAGE_CAP=50000 \
         METRIC_ACCUMULATION_LEVEL=ALL_READS;
         """



# check the jupyter notebook for downstream analyses
