
import glob
from datetime import datetime

configfile: "NGSanalyze.json"

#print(config["output_file_part2"]+'Normalized.xlsx')
#print(expand(config["translate_output_directory"]+"{name}.xlsx", name=NAMES))

translate_fastq_files = glob_wildcards(config["translate_input"]+'{name}.fastq')
NAMES = sorted(translate_fastq_files.name)

expand(config["translate_input"]+"{name}.fastq", name=NAMES)
now = datetime.now() # current date and time
date_time = now.strftime("%m%d%Y")


rule translate:
    input:
        config["translate_input"] if os.path.isfile(config["translate_input"]) else expand(config["translate_input"]+"{name}.fastq", name=NAMES)
    output:
        config["translate_output"] if os.path.isfile(config["translate_output"]) else expand(config["translate_output"]+str(date_time)+"{name}.xlsx", name=NAMES)
    script:
        "src/translatefastq2_12_24.py"

Pep_files = glob_wildcards(config["peptide_input"]+'{name}.xlsx')
pepNAMES = sorted(Pep_files.name)

rule peptide_qualities:
    input:
        config["peptide_input"] if os.path.isfile(config["peptide_input"]) else expand(config["peptide_input"]+"{name}.xlsx", name=NAMES)
    output:
        config["peptide_output"] if os.path.isfile(config["peptide_input"]) else expand(config["peptide_output"]+"{name}Characteristics.xlsx", name=NAMES)
    script:
        "src/PeptideQualities.py"

rule Fasta:
    input:
        config["Fasta_input"] if os.path.isfile(config["Fasta_input"]) else expand(config["Fasta_input"]+"{name}.xlsx", name=NAMES)
    output:
        config["Fasta_output"] if os.path.isfile(config["Fasta_input"]) else expand(config["Fasta_output"]+"{name}Characteristics.xlsx", name=NAMES)
    script:
        "src/MakeFasta.py"

rule AApercentage:
    input:
        config["Fasta"] if os.path.isfile(config["Fasta"]) else expand(config["Fasta"]+"{name}.fasta", name=NAMES)
    output:
        #config["Fasta_output"] if os.path.isfile(config["Fasta_input"]) else expand(config["Fasta_output"]+"{name}Characteristics.xlsx", name=NAMES)
    script:
        "src/AApercentage.py"

NGS_files = glob_wildcards(config["growth_directory"]+'{name}.xlsx')
NAMES = sorted(NGS_files.name)

rule growth:
    input:
        expand(config["growth_directory"]+"{name}.xlsx", name=NAMES)
    output:
        config["output_file_growth"]
    script:
        "src/Growth.py"

rule AddGrowth:
    input:
        expand(config["growth_file"])+expand(config["output_file_growth"])
    output:
        config["output_Score2"]
    script:
        "src/AddGrowth.py"

NGS_files = glob_wildcards(config["ForestInput"]+'{name}.xlsx')
NAMES = sorted(NGS_files.name)

rule Forest:
    input:
        expand(config["ForestInput"]+"{name}.xlsx", name=NAMES)
    output:
        config["ForestOutput"]
    script:
        "src/CNN.py"

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

#print(expand(config["pos_input_directory"]+"{name}.xlsx", name=NAMES))
rule design_choices:
    input:
        expand(config["pos_input_directory"]+"{name}.xlsx", name=NAMES) +expand(config["neg_input_directory"]+"{name2}.xlsx", name2=NAMES2)
    output:
        config["output_file_part2"]
    script:
        "src/NGSanalyzeMain3_26_24.py"

rule false_positives:
    input:
        expand(config["translate_output"]+"{name}.xlsx", name=NAMES)

rule visualization:
    input:
    script:
        "src/Visualize.py"
    #normalized score on y axis, (avg positive, and avg negative) diagonal not interesting

    #20in both positive libraries or 40 in 1 and 0 and the other
    #evenness of measurement
    #can rank by even-ness
        #if more even in positive, then more likely we care about it
            #care about max of negative library

rule motif_finding:
    input:
