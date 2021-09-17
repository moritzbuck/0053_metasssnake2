necessary_paths = ["raw_folder", "config_file", "libraries_file", "assemblies_file", "binnings_file", "binsets_file"]
optional_library_fields = ['fwd', 'rev', 'unp', 'read_subset', 'sourmash_k', 'sourmash_scaled', 'sortmerna_refs']

general_fields = {
    'seqio_buffer_size' : {
    'type' : int,
    'default' : 1000,
    'min' : 1,
    },
    'loglinelen' : {
    'type' : int,
    'default' : 90,
    'min' : 0
    },
    'other_parameters' : {
    'type' : str,
    'default' : ''
    },
    'temp_folder' : {
    'str' : str,
    'default' : '/tmp/'
    }
}

libraries_fields = {
    'fwd' : {
    'type' : str,
    'default' : ""
    },
    'rev' : {
    'type' : str,
    'default' : ""
    },
    'unp' : {
    'type' : str,
    'default' : ""
    },
    'read_subset' : {
    'type' : int,
    'default' : 100000,
    'min' : 0
    },
    'sourmash_k' : {
    'type' : int,
    'default' : 31,
    'min' : 11
    },
    'sourmash_scaled' : {
    'type' : int,
    'default' : 1000,
    'min' : 99
    },
    'sortmerna_refs' : {
    'type' : str,
    'default' : "~/dbs/sortmerna/set5-database.fasta;~/dbs/sortmerna/set6-database.fasta"
    },
    'other_parameters' : {
    'type' : str,
    'default' : {}
    }
}

necessary_assemblies_fields = ['libraries']
assemblies_fields =  {
    'libraries' : {
        'type' : str
    },
    'assembler' : {
        'possibles' : ['megahit', 'spades'],
        'default' : 'megahit'
    },
    'preprocess' : {
    'possibles' : [ 'bbnorm', 'none' ],
    'default' : 'none'
    },
    'length_cutoff' : {
    'type' : int,
    'default' : 2500,
    'min' : 0,
    },
    'other_parameters' : {
    'type' : str,
    'default' : ''
    },
    'keep_unpaired' : {
    'type' : bool,
    'default' : {}
    }
}

necessary_binnings_fields = ['libraries', 'assemblies']
binnings_fields =  {
    'libraries' : {
        'type' : str
    },
    'assemblies' : {
        'type' : str
    },
    'binner' : {
        'possibles' : ['metabat', 'vamb'],
        'default' : 'metabat'
    },
    'min_bin_size' : {
    'type' : int,
    'default' : 500000,
    'min' : 0,
    },
    'other_parameters' : {
    'type' : str,
    'default' : {}
    }
}

necessary_binsets_fields = ['binnings']
binsets_fields =  {
    'binnings' : {
        'type' : str
    },
    'other_parameters' : {
        'type' : str,
        'default' : {}
    },
    'min_completeness' : {
        'type' : float,
        'default' : 30,
        'min' : 0,
        'max' : 100
    },
    'max_contamination' : {
        'type' : float,
        'default' : 5,
        'min' : 0,
    },
    'min_coding' : {
        'type' : float,
        'default' : 0.2,
        'min' : 0,
        'max' : 1
    },
    'min_size' : {
        'type' : float,
        'default' : 200000,
        'min' : 0
    },
    'keep_fails' : {
        'type' : bool,
        'default' : True
    },
    'additional_annotation' : {
        'type' : str,
        'default' : "",
        'possibles' : ['gtdbtk+sourmash', 'eggnogmapper', 'anvioscg', 'motulize' ]
    },
}

necessary_mappings_fields = ['binset', 'libraries']
mappings_fields =  {
    'libraries' : {
        'type' : str
    },
    'binset' : {
        'type' : str
    },
    'mapper' : {
        'possibles' : ['bwa-mem2'],
        'default' : 'bwa-mem2'
    },
    'subset' : {
        'type' : bool,
        'default' : False
    },
    'taxfield' : {
    'type' : str,
    'default' : 'scg_taxo'
    },
    'other_parameters' : {
    'type' : str,
    'default' : {}
    },
    'min_nucleotide_id' : {
        'type' : float,
        'default' : 50,
        'min' : 0
    },
}
