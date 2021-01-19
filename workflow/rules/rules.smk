from os.path import join as pjoin

rule library_processing:
    input : unpack(lambda wildcards :  [pjoin(config['raw_folder'],v) for v in config['libraries'][wildcards.lib_name]['fwd'] + config['libraries'][wildcards.lib_name]['rev'] + config['libraries'][wildcards.lib_name]['unp'] ])
    output :
        fwd = "{root}/libraries/{lib_name}/{lib_name}_fwd.fastq.gz",
        rev = "{root}/libraries/{lib_name}/{lib_name}_rev.fastq.gz",
        unp = "{root}/libraries/{lib_name}/{lib_name}_unp.fastq.gz",
        sig = "{root}/libraries/{lib_name}/{lib_name}.sig.gz",
        sub_fwd = "{root}/libraries/{lib_name}/subs/subs_{lib_name}_fwd.fastq.gz",
        sub_rev = "{root}/libraries/{lib_name}/subs/subs_{lib_name}_rev.fastq.gz",
        sub_unp = "{root}/libraries/{lib_name}/subs/subs_{lib_name}_unp.fastq.gz",
        fastp_stats = "{root}/libraries/{lib_name}/stats/fastp_{lib_name}.stats",
    log :
        fastp_log = directory("{root}/libraries/{lib_name}/logs/fastp_logs/"),
        reformat_log = "{root}/libraries/{lib_name}/logs/bbtools_reformat.log",
        sourmash_log = "{root}/libraries/{lib_name}/logs/sourmash_compute.log",
        log = "{root}/libraries/{lib_name}/logs/{lib_name}.log",
        env = "{root}/libraries/{lib_name}/logs/library_processing.yaml",
        settings = "{root}/libraries/{lib_name}/settings/library_processing_settings.json"
    threads : 24
    params : script = "workflow/scripts/library_processing.py", config_file = pjoin(config['root_folder'],config['config_file'])
    conda : "../envs/library_processing.yaml"
    shell : """
        python {params.script} {wildcards.lib_name} {params.config_file} {wildcards.root} {wildcards.root}/libraries/{wildcards.lib_name}/ {threads}
        """

rule library_rrna_spliting:
    input :
        fwd = "{root}/libraries/{lib_name}/{lib_name}_fwd.fastq.gz",
        rev = "{root}/libraries/{lib_name}/{lib_name}_rev.fastq.gz",
        unp = "{root}/libraries/{lib_name}/{lib_name}_unp.fastq.gz",
    output :
        fwd_mrna = "{root}/libraries/{lib_name}/{lib_name}_mrna_fwd.fastq.gz",
        rev_mrna = "{root}/libraries/{lib_name}/{lib_name}_mrna_rev.fastq.gz",
        unp_mrna = "{root}/libraries/{lib_name}/{lib_name}_mrna_unp.fastq.gz",
        fwd_rrna = "{root}/libraries/{lib_name}/{lib_name}_rrna_fwd.fastq.gz",
        rev_rrna = "{root}/libraries/{lib_name}/{lib_name}_rrna_rev.fastq.gz",
        unp_rrna = "{root}/libraries/{lib_name}/{lib_name}_rrna_unp.fastq.gz",
        sortmerna_stats_paired = "{root}/libraries/{lib_name}/stats/sortmerna_{lib_name}_paired.stats",
        sortmerna_stats_unpaired = "{root}/libraries/{lib_name}/stats/sortmerna_{lib_name}_unpaired.stats"
    log :
        sortmerna_log = "{root}/libraries/{lib_name}/logs/sortmerna.log",
        log = "{root}/libraries/{lib_name}/logs/{lib_name}.log",
        env = "{root}/libraries/{lib_name}/logs/library_rrna_spliting.yaml",
        settings = "{root}/libraries/{lib_name}/logs/library_rrna_spliting_settings.json"
    threads : 24
    params : script = "workflow/scripts/sortmerna_wrapper.py", config_file = pjoin(config['root_folder'],config['config_file'])
    conda : "../envs/sortmerna.yaml"
    shell : """
        python {params.script} {wildcards.lib_name} {params.config_file} {wildcards.root} {wildcards.root}/libraries/{wildcards.lib_name}/ {threads}
        """

rule assembly:
    input : unpack(lambda wildcards :  [ "{root}/libraries/{lib_name}/{lib_name}_fwd.fastq.gz".format(root = wildcards.root, lib_name = lib) for lib in config['assemblies'][wildcards.ass_name]['libraries'] ])
    output :
        assembly = "{root}/assemblies/{ass_name}/assembly.fna"
    log :
        megahit_log = "{root}/assemblies/{ass_name}/logs/megahit.log",
        log = "{root}/assemblies/{ass_name}/logs/{ass_name}.log",
        env = "{root}/assemblies/{ass_name}/logs/assembly.yaml",
        settings = "{root}/assemblies/{ass_name}/logs/assembly_settings.json"
    threads : 24
    params : script = "workflow/scripts/assemble.py", config_file = pjoin(config['root_folder'],config['config_file'])
    conda : "../envs/assembly.yaml"
    shell : """
        python {params.script} {wildcards.ass_name} {params.config_file} {wildcards.root} {wildcards.root}/assemblies/{wildcards.ass_name}/ {threads}
        """

rule binning:
    input : unpack(lambda wildcards :  [ "{root}/libraries/{lib_name}/{lib_name}_fwd.fastq.gz".format(root = wildcards.root, lib_name = lib) for lib in config['binnings'][wildcards.binning_name]['libraries'] ]),
            unpack(lambda wildcards :  [ "{root}/assemblies/{ass}/assembly.fna".format(root = wildcards.root, ass = ass) for ass in config['binnings'][wildcards.binning_name]['assemblies'] ])
    output :
        assembly_by_bins = "{root}/binnings/{binning_name}/binned_assembly.fna"
    log :
        log = "{root}/binnings/{binning_name}/logs/binning.log",
        env = "{root}/binnings/{binning_name}/logs/binning.yaml",
        settings = "{root}/binnings/{binning_name}/logs/binning_settings.json"
    threads : 24
    params : script = "workflow/scripts/binning.py", config_file = config['config_file']
    conda : "../envs/binning.yaml"
    shell : """
        python {params.script} {wildcards.binning_name} {params.config_file} {wildcards.root} {wildcards.root}/binnings/{wildcards.binning_name}/ {threads}
        """

rule binsetting:
    input : unpack(lambda wildcards :  [ "{root}/binnings/{binning}/binned_assembly.fna".format(root = wildcards.root, binning = binning) for binning in config['binsets'][wildcards.binset_name]['binnings'] ])
    output :
        assembly_by_binset = "{root}/binsets/{binset_name}/binset.fna"
    log :
        log = "{root}/binset/{binset_name}/logs/binset.log",
        env = "{root}/binnings/{binset_name}/logs/binset.yaml",
        settings = "{root}/binnings/{binset_name}/logs/binset_settings.json"
    threads : 24
    params : script = "workflow/scripts/binset.py", config_file = config['config_file']
    conda : "../envs/binset.yaml"
    shell : """
        python {params.script} {wildcards.binset_name} {params.config_file} {wildcards.root} {wildcards.root}/binnings/{wildcards.binset_name}/ {threads}
        """
