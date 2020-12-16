assembly_fields =  {
    'libraries' : {
        'type' : list
    },
    'bin_mapping' : {
            'type' : list,
            'default' : []
    },
    'assembler' : {
        'possibles' : ['megahit'],
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
    'seqio_buffer' : {
    'type' : int,
    'default' : 1000,
    'min' : 0
    }
}

libraries_fields = {
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
    'type' : list,
    'default' : ["~/dbs/sortmerna/set5-database.fasta", "~/dbs/sortmerna/set6-database.fasta"]
    }
}

rna_libraries = {
'default' : [],
'type' : list
}



loglinelen = 90
