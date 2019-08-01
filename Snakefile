
# "rules" at the top are basically all the important files that you want output when 
# snakemake is run, these are the "requirements" that snakemake looks for when it is run

fastqc_output = ["data/0Hour_001_1_fastqc.html", "data/6Hour_001_1_fastqc.html",
  "data/0Hour_001_2_fastqc.html", "data/6Hour_001_2_fastqc.html",
  "data/0Hour_002_1_fastqc.html", "data/6Hour_002_1_fastqc.html",
  "data/0Hour_002_2_fastqc.html", "data/6Hour_002_2_fastqc.html"]

# A common error is to forget the comma between the input or output items. Since Python
# concatenates subsequent strings, this can lead to unexpected behavior.

rule all:
  input:
    "multiqc_report.html"
    
rule fastqc_a_file:
  input:
    "{filename}.fq.gz"
  output:
    "{filename}_fastqc.html",
    "{filename}_fastqc.zip"
  shell:
    "fastqc {input}"
    
rule run_multiqc:
  input:
    fastqc_output
  output:
    "multiqc_report.html",
    directory("multiqc_data")
  shell:
    "multiqc data/"

#rule trim_reads:
#  input:
#    "{filename}_1.fq.gz",
#    "{filename}_2.fq.gz"
#  output:
#    "{filename}_1.pe.qc.fq.gz",
#    "{filename}_1.se.qc.fq.gz",
#    "{filename}_2.pe.qc.fq.gz",
#    "{filename}_2.se.qc.fq.gz"
#  shell:
#    """trimmomatic PE {input} {output} LEADING:2 TRAILING:2 \
#      SLIDINGWINDOW:4:15 \
#      MINLEN:25"""

# will automatically delete stuff for us

#rule clean:
#  shell:
#    "rm -f {fastqc_output} multiqc_report.html"

# more rules for later

# rule fastqc_a_file

# rule trimm_a_file

# rule hisat2_a_file

# rule featurecount_a_file

# rule tximeta_a_file 
