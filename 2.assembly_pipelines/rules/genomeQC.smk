def get_meryl_union_sum_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Specimen"] == wildcards.specimen]
    samples = samples.to_records(index=False)
    return ["output/meryl/{specimen}/{settings}/{facility}/{lane}/{id}.meryl".format(specimen=s[0],
            settings = wildcards.settings, facility = s[1], lane = s[2], id = s[3]) for s in samples]


def get_jellyfish_inputs(wildcards):
    samples = pd.read_table("pepsamples.tsv", index_col=False)
    samples = samples[samples["Specimen"] == wildcards.specimen]
    samples = samples.to_records(index=False)
    returnlist = ["output/HiFi-adapterFiltered/{specimen}/{settings}/{facility}/{lane}/{id}.ccs.filt.fastq.gz".format(specimen=s[0],
                  settings = wildcards.settings, facility = s[1], lane = s[2], id = s[3]) for s in samples]
    #print(returnlist)
    return(returnlist)


rule meryl_run:
    input: "output/HiFi-adapterFiltered/{specimen}/{settings}/{facility}/{lane}/{id}.ccs.filt.fastq.gz"
    output: directory("output/meryl/{specimen}/{settings}/{facility}/{lane}/{id}.meryl")
    log: "logs/meryl/{specimen}/{settings}/{facility}/{lane}/{id}.log"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "meryl --version && meryl count threads={threads} memory=300 k=21 {input} output {output} &> {log}"

rule meryl_union_sum:
    input: get_meryl_union_sum_inputs
    output: directory("output/meryl/union-sum/{specimen}/{settings}/{specimen}.meryl")
    log: "logs/meryl-union-sum/{specimen}/{settings}/{specimen}.log"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "meryl union-sum threads={threads} {input} output {output} &> {log}"
    
rule merqury_run_haplotig:
    input: 
        meryl = "output/meryl/union-sum/{specimen}/{settings}/{specimen}.meryl",
        hap1_p = "output/hifiasm/{specimen}/{settings}/{opt}/{specimen}.asm.bp.hap1.p_ctg.gfa",
        hap2_p = "output/hifiasm/{specimen}/{settings}/{opt}/{specimen}.asm.bp.hap2.p_ctg.gfa"        
    output: directory("output/merqury/{specimen}/{settings}/{opt}/{specimen}.haplotig")
    conda: "../envs/genomeQC.yml"
    threads: 52
    # shadow: "shallow"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        echo $ROOTPROJDIR
        mkdir -p {output}
        cd {output}
        pwd -P
        #rmdir $ROOTPROJDIR/{output}
        merqury.sh $ROOTPROJDIR/{input.meryl} $ROOTPROJDIR/{input.hap1_p} $ROOTPROJDIR/{input.hap2_p} {wildcards.specimen}.haplotig
        cd -
    """

rule merqury_run_consensus:
    input: 
        meryl = "output/meryl/union-sum/{specimen}/{settings}/{specimen}.meryl",
        consensus = "output/hifiasm/{specimen}/{settings}/{opt}/{specimen}.asm.bp.p_ctg.gfa",
    output: directory("output/merqury/{specimen}/{settings}/{opt}/{specimen}.consensus")
    conda: "../envs/genomeQC.yml"
    threads: 52
    # shadow: "shallow"
    shell: """
        ROOTPROJDIR="$(pwd -P)"
        echo $ROOTPROJDIR
        mkdir -p {output}
        cd {output}
        pwd -P
        #rmdir $ROOTPROJDIR/{output}
        merqury.sh $ROOTPROJDIR/{input.meryl} $ROOTPROJDIR/{input.consensus} {wildcards.specimen}.consensus
        cd -
    """


rule jellyfish_run:
    input: get_jellyfish_inputs
    output: 
        bc = "output/jellyfish/{specimen}/{settings}/{specimen}.bc",
        jf = "output/jellyfish/{specimen}/{settings}/{specimen}.jf"
    conda: "../envs/genomeQC.yml"
    version: "0.0.2"
    threads: 52
    params:
        kmer = "21",
        generatorfile = lambda wildcards: "generators_{a}_{b}".format(a=wildcards.specimen, b=wildcards.settings)
    shell: """
        ls {input} | xargs -I[] echo zcat [] > {params.generatorfile}
        jellyfish bc -m {params.kmer} -s 100G -g {params.generatorfile} -G {threads} -t {threads} -o {output.bc} 
        jellyfish count -C -m {params.kmer} -s 3G -g {params.generatorfile} -G {threads} -t {threads} -o {output.jf} 
        rm {params.generatorfile}
    """

rule jellyfish_histo:
    input: "output/jellyfish/{specimen}/{settings}/{specimen}.jf"
    output: "output/jellyfish/{specimen}/{settings}/{specimen}.histo"
    conda: "../envs/genomeQC.yml"
    threads: 52
    shell: "jellyfish histo -t {threads} {input} > {output}"


rule genomescope:
    input: "output/jellyfish/{specimen}/{settings}/{specimen}.histo"
    output: 
        linearplot = "output/genomescope2/{specimen}/{settings}/{specimen}/linear_plot.png",
        logplot = "output/genomescope2/{specimen}/{settings}/{specimen}/log_plot.png",
        model = "output/genomescope2/{specimen}/{settings}/{specimen}/model.txt",
        progress = "output/genomescope2/{specimen}/{settings}/{specimen}/progress.txt",
        summary = "output/genomescope2/{specimen}/{settings}/{specimen}/summary.txt",
        translinearplot = "output/genomescope2/{specimen}/{settings}/{specimen}/transformed_linear_plot.png",
        translogplot = "output/genomescope2/{specimen}/{settings}/{specimen}/transformed_log_plot.png"
    params:
       prefix = "output/genomescope2/{specimen}/{settings}/{specimen}",
       kmer = 21
    log: "output/genomescope2/{specimen}/{settings}/{specimen}/log.txt"
    conda: "../envs/genomeQC.yml"
    shell: "which R && genomescope2 -i {input} -o {params.prefix} -k {params.kmer} > {log} 2>&1"

def get_genomesize(gs2summary):
    import re
    genomesize = 0 
    pat = re.compile("(?<=len:)\d+")
    with open(gs2summary) as infile: 
        for line in infile:
            if pat.findall(line):
                genomesize = pat.search(line).group(0) 
    if genomesize: 
        return genomesize
    else:
        raise Exception("no genomesize found - did the file output from genomescope2 change?")

rule genometools:
    input: 
        fasta = "output/hifiasm-fasta/{specimen}/{settings}/{opt}/{specimen}.{genometype}.fa",
        gs2summary = "output/genomescope2/{specimen}/{settings}/{specimen}/log.txt"
    #output: "output/gt-seqstat/{specimen}/{settings}/{opt}/{specimen}.{genometype}.stats"
    output: "output/gt-seqstat/{settings}/{opt}/{specimen}.{genometype}.stats"
    run: 
        shell("./code/bin/gt "
              "seqstat -contigs -genome "
              "{genomesize} "
              "{input} > "
              "{output}".format(input=input.fasta,
                                output=output, 
                                genomesize=get_genomesize(input.gs2summary))
              )

    
