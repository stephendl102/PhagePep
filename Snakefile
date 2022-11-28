import glob

configfile: "NGSanalyze.json"




#print(config["output_file_part2"]+'Normalized.xlsx')
#print(expand(config["translate_output_directory"]+"{name}.xlsx", name=NAMES))

translate_fastq_files = glob_wildcards(config["translate_input_directory"]+'{name}.fastq')
NAMES = sorted(translate_fastq_files.name)

rule translate_fastq:
    input:
        config["translate_input_directory"] if os.path.isfile(config["translate_input_directory"]) else expand(config["translate_input_directory"]+"{name}.fastq", name=NAMES)
    output:
        config["translate_output_directory"] if os.path.isfile(config["translate_input_directory"]) else expand(config["translate_output_directory"]+"{name}.fastq", name=NAMES)
    script:
        "src/translatefastq2.py"

Pep_files = glob_wildcards(config["peptide_input"]+'{name}.xlsx')
pepNAMES = sorted(Pep_files.name)

rule peptide_qualities:
    input:
        config["peptide_input"] if os.path.isfile(config["peptide_input"]) else expand(config["peptide_input"]+"{name}.xlsx", name=NAMES)
    output:
        config["peptide_output"] if os.path.isfile(config["peptide_input"]) else expand(config["peptide_output"]+"{name}.fa", name=NAMES)
    script:
        "src/PeptideQualities.py"

rule sql:
    input:
        config["SQL_input_directory"] if os.path.isfile(config["SQL_input_directory"]) else expand(config["SQL_input_directory"]+"{name}.xlsx", name=NAMES)
    output:
        expand(config["SQL_output_directory"]+"Database.sql")

NGS_files = glob_wildcards(config["pos_input_directory"]+'{name}.xlsx')
NAMES = sorted(NGS_files.name)

NGS_files2 = glob_wildcards(config["neg_input_directory"]+'{name2}.xlsx')
NAMES2 = sorted(NGS_files2.name2)
#print(NAMES2)


rule design_choices:
    input:
        expand(config["pos_input_directory"]+"{name}.xlsx", name=NAMES)+expand(config["neg_input_directory"]+"{name2}.xlsx", name2=NAMES2)
    output:
        expand(config["output_file_part2"])
    script:
        "src/NGSanalyzeMain2.py"

rule false_positives:
    input:
        expand(config["translate_output_directory"]+"{name}.xlsx", name=NAMES)

rule visualization:
    input:
        expand(config["output_file_part2"]+"Normalized.xlsx")
    #normalized score on y axis, (avg positive, and avg negative) diagonal not interesting

    #20in both positive libraries or 40 in 1 and 0 and the other
    #evenness of measurement
    #can rank by even-ness
        #if more even in positive, then more likely we care about it
            #care about max of negative library

rule motif_finding:
    input:
