rule bbduk:
    input:
        sample = [dir_in_trim + "/{sample}_{seq_lane}.fastq.gz"]
    output:
        trimmed = "results/trimmed/{sample}/{sample}_{seq_lane}_trimmed.fastq.gz",
        singleton = "results/trimmed/{sample}/{sample}_{seq_lane}_single.fastq.gz",
        discarded = "results/trimmed/{sample}/{sample}_{seq_lane}_discarded.fastq.gz",
        stats = "results/trimmed/{sample}/{sample}_{seq_lane}_stats.txt"
    conda: 
        config["conda_envs"]["preprocessing"]
    threads: 
        get_resource(config, "bbduk", "threads")
    resources:
        mem_mb = get_resource(config, "bbduk", "mem_mb"),
        runtime = get_resource(config, "bbduk", "runtime")
    params:
        adapters = config["adapters"],
        polyA = config["polyA"],
        extra = config["parameters"]["bbudk"]["extra"]
    log:
        "log/bbduk/{sample}_{seq_lane}_bbduk.log"
    benchmark:
        "benchmarks/bbduk/{sample}_{seq_lane}_bbduk.bmk"
    shell:
        """
        bbduk.sh in={input.sample} \
            out={output.trimmed} \
            outs={output.singleton} \
            outm={output.discarded} \
            stats={output.stats} \
            ref={params.adapters},{params.polyA} \
            threads={threads} \
            {params.extra} > {log} 2>&1
        """
