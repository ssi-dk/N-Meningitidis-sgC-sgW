rule all:
    input:
        auspice_json = "auspice/mg.json",

# Config variables to be used by rules
# Parameters are defined within their own rules

rule config:
    params:
        seq = "/Volumes/ssDisk/Projects/1908_SgW_meningokokker/bestsnp_matrices/bestsnp_N2.filtered_polymorphic_sites_manal_purge_N-_stripped_wN.fasta",
        ref = "in/reference.fasta",
        meta = "in/metadata.tsv",
        exclude = "config/dropped_strains.txt",
        cores = "8",
        colors = "config/color.tsv",
        geo_info = "config/my_lat_long.tsv",
        weights = 'config/weights.tsv',
        config = "config/config.json"

config = rules.config.params #so we can use config.x rather than rules.config.params.x
#end of config definition

rule filter:
    input:
        seq = config.seq,
        meta = config.meta,
        exclude = config.exclude
    output:
        "results/filtered.fasta"
    shell:
        """
        augur filter --sequences {input.seq} \
            --metadata {input.meta} \
            --exclude {input.exclude} \
            --output {output}
        """

rule tree:
    input:
        aln = rules.filter.output,
        ref = config.ref,
        # sites = config.sites
    output:
        "results/tree_raw.nwk"
    params:
        method = 'iqtree',
        model = 'GTR', # According to ModelTest in independent IQTREE run
        cores = config.cores
    shell:
        """
        augur tree --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --substitution-model {params.model} \
            --nthreads {params.cores} \
            --output {output}
        """

rule refine:
    input:
        tree = rules.tree.output,
        aln = rules.filter.output,
        metadata = config.meta,
        ref = config.ref
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json",
    params:
        root = 'best',
        coal = 'opt',
        # clock_rate = "0.000220",
        # clock_sd = "0.000291",
        date_inference = 'joint', 
        br_inference = 'input'
    shell:
            # --clock-rate {params.clock_rate} \
            # --clock-std-dev {params.clock_sd} \
        """
        augur refine --tree {input.tree} \
            --alignment {input.aln} \
            --metadata {input.metadata} \
            --coalescent {params.coal} \
            --root {params.root} \
            --timetree \
            --date-confidence \
            --keep-polytomies \
            --date-inference {params.date_inference} \
            --branch-length-inference {params.br_inference} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.filter.output,
        ref = config.ref
    output:
        nt_data = "results/nt_muts.json",
        # vcf_out = "results/nt_muts.vcf"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} \
            --alignment {input.alignment} \
            --vcf-reference {input.ref} \
            --keep-ambiguous \
            --inference {params.inference} \
            --output-node-data {output.nt_data}
        """

rule traits:
    input:
        # tree = "results/tree.nwk",
        tree = rules.refine.output.tree,
        meta = config.meta,
        weights = config.weights
    output:
        "results/traits.json"
        # rules.all.input.tmp_out
    params:
        traits = 'serogroup country continent cluster 2013_strain Lineages_Eriksson_2018 Lineages_Lucidarme_2015 Lineages_Knol_2017'
    shell:
        """
        augur traits --tree {input.tree} \
            --metadata {input.meta} \
            --weights {input.weights} \
            --columns {params.traits} \
            --output-node-data {output}
        """

rule export:
    input:
        # tree = "results/tree.nwk",
        # branch_lengths = "results/branch_lengths.json",
        # traits = "results/traits.json",
        # nt_muts = "results/nt_muts.json",
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        metadata = config.meta,
        color_defs = config.colors,
        config = config.config,
        geo_info = config.geo_info,
    output:
        out_json = rules.all.input.auspice_json
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} \
            --auspice-config {input.config} \
            --colors {input.color_defs} \
            --lat-longs {input.geo_info} \
            --output {output.out_json}
        # augur validate --json {output.out_json}
        """
