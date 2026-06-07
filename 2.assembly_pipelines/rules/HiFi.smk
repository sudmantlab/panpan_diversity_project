import os
import pandas as pd

rule HiFiAdapterFilt:
    input: "PANPAN_ccs/{specimen}/{facility}/{lane}/{id}.ccs.bam" 
    output: 
        filtered = "output/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.ccs.filt.fastq.gz",
        stats = "output/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.ccs.stats",
    log: "logs/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.log"
    version: 2.00
    params:
        outDir = lambda wildcards, output: os.path.dirname(output[0]),
        inBaseName = lambda wildcards, input: os.path.basename(input[0]),
        inDir = lambda wildcards, input: os.path.dirname(input[0]),
        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0],        
    conda: "../envs/HiFiAssembly.yml"
    threads: 10
    shell: """
        ROOTPROJDIR=$(pwd -P)
        cd {params.inDir}
        #mkdir $ROOTPROJDIR/{params.outDir}
        $ROOTPROJDIR/code/HiFiAdapterFilt/pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
        cd -
    """

#rule HiFiAdapterFilt_fastq:
#    input: "PANPAN_ccs/{specimen}/{facility}/{lane}/{id}.ccs.fastq.gz" 
#    output: 
#        filtered = "output/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.ccs.filt.fastq.gz",
#        stats = "output/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.ccs.stats",
#    log: "logs/HiFi-adapterFiltered/{specimen}/joana_settings_shortcut/{facility}/{lane}/{id}.log"
#    #version: 2.00
#    params:
#        outDir = lambda wildcards, output: os.path.dirname(output[0]),
#        inBaseName = lambda wildcards, input: os.path.basename(input[0]),
#        inDir = lambda wildcards, input: os.path.dirname(input[0]),
#        inPref = lambda wildcards, input: os.path.splitext(os.path.basename(input[0]))[0].replace(".fastq","")
#    conda: "../envs/HiFiAssembly.yml"
#    threads: 10
#    shell: """
#        ROOTPROJDIR=$(pwd -P)
#        echo $ROOTPROJDIR
#        cd {params.inDir}
#        mkdir -p $ROOTPROJDIR/{params.outDir}
#        $ROOTPROJDIR/code/HiFiAdapterFilt/pbadapterfilt.sh -p {params.inPref} -t {threads} -o $ROOTPROJDIR/{params.outDir} &> $ROOTPROJDIR/{log}
#        cd -
#        """


def get_hifiasm_inputs(wildcards):
    hifi_path = "output/HiFi-adapterFiltered/{specimen}/{settings}/{facility}/{lane}/{id}.ccs.filt.fastq.gz"
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Specimen"] == wildcards.specimen]
    samples = samples.to_records(index=False)
    input_samples = [hifi_path.format(specimen=s[0], settings = wildcards.settings, facility = s[1], lane = s[2], id = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check pepsamples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples



rule hifiasm_noopts:
    version: "0.2.2"
    input: get_hifiasm_inputs
    output: 
        p_ctg_hap1 = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_ctg.gfa",
        p_ctg_lowQ = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_ctg.noseq.gfa",
        p_utg = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_utg.gfa",
        p_utg_lowQ = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_utg.lowQ.bed",
        p_utg_noseq = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.p_utg.noseq.gfa",
        r_utg = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.r_utg.gfa",
        r_utg_lowQ = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.r_utg.lowQ.bed",
        r_utg_noseq = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.bp.r_utg.noseq.gfa",
        ec = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.ec.bin",
        ovlp_reverse = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm"
    threads: 52
    log: "output/hifiasm/{specimen}/{settings}/no_opts/{specimen}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    # shell: "hifiasm -l2 -o {params.prefix} -t {threads} {input} > {log} 2>&1"
    shell: "hifiasm -o {params.prefix} -t {threads} {input} > {log} 2>&1" 


def get_hifiasm_hic_inputs(wildcards):
    file_path = path_trimmed = "output/trimmed-hic/{specimen}/{sample_name}-trimmed_{read}.fastq.gz"
    input_dict = {"left": [], "right": [], "hifi": get_hifiasm_inputs(wildcards)}
    samples = pd.read_table("pepsamples_hic.tsv", sep= '\t',index_col=False)
    samples = samples[samples.specimen==wildcards.specimen][samples.type == "Hi-C"]
    print(samples)
    samples_grouped = samples.groupby(samples.sample_name)
    for sample in set(samples["sample_name"].tolist()):
        sample_subset = samples_grouped.get_group(sample)
        if len(sample_subset) == 0:
           raise Exception("No files available for sample {}".format(sample))
        input_dict["left"].extend([path_trimmed.format(specimen = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R1"].itertuples()])
        input_dict["right"].extend([path_trimmed.format(specimen = r[2], sample_name = r[3], read = r[4]) for r in sample_subset[sample_subset["read"] == "R2"].itertuples()])
    return input_dict


rule hifiasm_noopts_HiC:
    version: "0.2.1"
    input: unpack(get_hifiasm_hic_inputs)
    output:
        p_ctg_hap1 = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap1.p_ctg.gfa",
        p_ctg_hap1_lowQ = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap1.p_ctg.lowQ.bed",
        p_ctg_hap1_noseq = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap1.p_ctg.noseq.gfa",
        p_ctg_hap2 = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap2.p_ctg.gfa",
        p_ctg_hap2_lowQ = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap2.p_ctg.lowQ.bed",
        p_ctg_hap2_noseq = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.hap2.p_ctg.noseq.gfa",
        p_ctg = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_ctg.gfa",
        p_ctg_lowQ = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_ctg.lowQ.bed",
        p_ctg_noseq = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_ctg.noseq.gfa",
        p_utg = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_utg.gfa",
        p_utg_lowQ = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_utg.lowQ.bed",
        p_utg_noseq = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.p_utg.noseq.gfa",
        r_utg = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.r_utg.gfa",
        r_utg_lowQ = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.r_utg.lowQ.bed",
        r_utg_noseq = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.r_utg.noseq.gfa",
        ec = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.ec.bin",
        #lk_bin = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.hic.lk.bin",
        ovlp_reverse = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.ovlp.reverse.bin",
        ovlp_source = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.ovlp.source.bin"
    params:
        prefix = "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm"
    threads: 32
    log: "output/hifiasm-HiC/{specimen}/{settings}/no_opts/{specimen}.asm.log"
    conda: "../envs/HiFiAssembly.yml"
    # shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} -l2 {input.hifi} > {log} 2>&1"
    shell: "hifiasm -o {params.prefix} --h1 {input.left} --h2 {input.right} -t {threads} {input.hifi} > {log} 2>&1"


rule gfaToFa:
    input: "output/hifiasm/{specimen}/{settings}/{opt}/{specimen}.asm.bp.{genometype}.gfa"
    output: "output/hifiasm-fasta/{specimen}/{settings}/{opt}/{specimen}.{genometype}.fa"
    log: "logs/gfaToFa/{specimen}/{settings}/{opt}/{specimen}.{genometype}.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "gfatools gfa2fa {input} > {output} 2> {log}"

rule gfaToFa_hic:
    input: "output/hifiasm-HiC/{specimen}/{settings}/{opt}/{specimen}.asm.hic.{genometype}.gfa"
    output: "output/hifiasm-fasta-HiC/{specimen}/{settings}/{opt}/{specimen}.{genometype}.hic.fa",
    log: "logs/gfaToFa-HiC/{specimen}/{settings}/{opt}/{specimen}.{genometype}.log"
    conda: "../envs/HiFiAssembly.yml"
    shell: "gfatools gfa2fa {input} > {output} 2> {log}"

