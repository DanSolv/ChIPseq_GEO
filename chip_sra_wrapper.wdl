version 1.0

# ─────────────────────────────────────────────────────────────────────────────
# chip_sra_wrapper.wdl
#
# Wrapper around the ENCODE ChIP-seq pipeline (chip.wdl) that downloads
# samples from SRA before running the analysis. Built for small batches:
# 1 input control + 1–4 IP replicates.
#
# Usage:
#   1. Fill in chip_sra_wrapper_inputs.json with your SRR accessions.
#   2. Run on Terra alongside chip.wdl (imported below).
#   3. Outputs land in the same place as a normal chip.wdl run.
#
# The wrapper only handles download + handoff. All pipeline parameters
# (aligner, genome_tsv, peak caller, etc.) are passed through unchanged.
# ─────────────────────────────────────────────────────────────────────────────

import "https://raw.githubusercontent.com/ENCODE-DCC/chip-seq-pipeline2/refs/heads/master/chip.wdl" as chip

struct SampleInfo {
    String accession
    Boolean paired_end
    Boolean is_control
    String label
}

workflow chip_sra_wrapper {

    meta {
        author: "DanSolv"
        email: ""
        description: "Downloads FASTQs from SRA and runs the ENCODE ChIP-seq pipeline. Supports 1 input control and 1-4 IP replicates."
    }

    input {
        Array[SampleInfo] samples

        File?   genome_tsv
        String  pipeline_type
        String  title       = "Untitled"
        String  description = "No description"

        String  aligner       = "bowtie2"
        Boolean align_only    = false
        Boolean true_rep_only = false
        Float   pval_thresh   = 0.01
        Float   idr_thresh    = 0.05

        String  docker      = "encodedcc/chip-seq-pipeline:v2.2.2"
        String  singularity = "https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/chip-seq-pipeline_v2.2.2.sif"
        String  conda       = "encd-chip"

        Int download_cpu     = 6
        Int download_disk_gb = 120
        Int download_mem_gb  = 8
    }

    # ── Download all samples in parallel ────────────────────────────────────
    scatter (s in samples) {
        call sra_to_fastq {
            input:
                accession   = s.accession,
                paired_end  = s.paired_end,
                cpu         = download_cpu,
                mem_gb      = download_mem_gb,
                disk_gb     = download_disk_gb,
                docker      = docker,
                singularity = singularity,
                conda       = conda
        }
    }

    # ── Split downloaded FASTQs into IP replicates vs control ───────────────
    call split_samples {
        input:
            accessions  = samples,
            r1_files    = sra_to_fastq.fastq_R1,
            r2_files    = sra_to_fastq.fastq_R2,
            docker      = docker,
            singularity = singularity,
            conda       = conda
    }

    # ── Hand off to the ENCODE chip pipeline ────────────────────────────────
    call chip.chip {
        input:
            title               = title,
            description         = description,
            pipeline_type       = pipeline_type,
            genome_tsv          = genome_tsv,

            paired_end          = split_samples.ip_paired_end,
            ctl_paired_end      = split_samples.ctl_paired_end,

            fastqs_rep1_R1      = if length(split_samples.ip_r1) > 0
                                    then [split_samples.ip_r1[0]] else [],
            fastqs_rep1_R2      = if length(split_samples.ip_r2) > 0
                                    then [split_samples.ip_r2[0]] else [],

            fastqs_rep2_R1      = if length(split_samples.ip_r1) > 1
                                    then [split_samples.ip_r1[1]] else [],
            fastqs_rep2_R2      = if length(split_samples.ip_r2) > 1
                                    then [split_samples.ip_r2[1]] else [],

            fastqs_rep3_R1      = if length(split_samples.ip_r1) > 2
                                    then [split_samples.ip_r1[2]] else [],
            fastqs_rep3_R2      = if length(split_samples.ip_r2) > 2
                                    then [split_samples.ip_r2[2]] else [],

            fastqs_rep4_R1      = if length(split_samples.ip_r1) > 3
                                    then [split_samples.ip_r1[3]] else [],
            fastqs_rep4_R2      = if length(split_samples.ip_r2) > 3
                                    then [split_samples.ip_r2[3]] else [],

            ctl_fastqs_rep1_R1  = if length(split_samples.ctl_r1) > 0
                                    then [split_samples.ctl_r1[0]] else [],
            ctl_fastqs_rep1_R2  = if length(split_samples.ctl_r2) > 0
                                    then [split_samples.ctl_r2[0]] else [],

            aligner             = aligner,
            align_only          = align_only,
            true_rep_only       = true_rep_only,
            pval_thresh         = pval_thresh,
            idr_thresh          = idr_thresh,

            docker              = docker,
            singularity         = singularity,
            conda               = conda
    }

    output {
        File    qc_report         = chip.report
        File    qc_json           = chip.qc_json
        Boolean qc_json_ref_match = chip.qc_json_ref_match
        Array[File] fastqs_R1     = sra_to_fastq.fastq_R1
        Array[File] fastqs_R2     = sra_to_fastq.fastq_R2
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: download one SRA accession and convert to FASTQ.gz
# ─────────────────────────────────────────────────────────────────────────────

task sra_to_fastq {
    input {
        String  accession
        Boolean paired_end
        Int     cpu
        Int     mem_gb
        Int     disk_gb
        String  docker
        String  singularity
        String  conda
    }

    command <<<
        set -euo pipefail

        echo "Downloading ~{accession} ..."

        prefetch \
            --max-size 50GB \
            --output-directory . \
            ~{accession}

        echo "Converting to FASTQ ..."
        fasterq-dump \
            ~{accession} \
            --split-files \
            --skip-technical \
            --threads ~{cpu} \
            --outdir .

        gzip -f *.fastq

        # Rename to predictable names
        if [ -f ~{accession}_1.fastq.gz ]; then
            mv ~{accession}_1.fastq.gz ~{accession}_R1.fastq.gz
        elif [ -f ~{accession}.fastq.gz ]; then
            mv ~{accession}.fastq.gz ~{accession}_R1.fastq.gz
        fi

        if [ -f ~{accession}_2.fastq.gz ]; then
            mv ~{accession}_2.fastq.gz ~{accession}_R2.fastq.gz
        fi

        # Always create R2 so the output declaration is satisfied.
        # For SE samples this will be an empty file and is ignored downstream.
        if [ ! -f ~{accession}_R2.fastq.gz ]; then
            touch ~{accession}_R2.fastq.gz
        fi

        echo "Done."
        ls -lh *.fastq.gz
    >>>

    output {
        File fastq_R1 = "~{accession}_R1.fastq.gz"
        File fastq_R2 = "~{accession}_R2.fastq.gz"
    }

    runtime {
        cpu:         cpu
        memory:      "~{mem_gb} GB"
        disks:       "local-disk ~{disk_gb} SSD"
        docker:      "ncbi/sra-tools:3.0.0"
        singularity: singularity
        conda:       "sra-tools"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: split scatter outputs into IP vs control arrays
# ─────────────────────────────────────────────────────────────────────────────

task split_samples {
    input {
        Array[SampleInfo] accessions
        Array[File]       r1_files
        Array[File]       r2_files
        String            docker
        String            singularity
        String            conda
    }

    File accessions_json = write_json(accessions)

    command <<<
        python3 <<'PYEOF'
import json

with open("~{accessions_json}") as f:
    samples = json.load(f)

r1_paths = "~{sep='|' r1_files}".split("|")
r2_paths = "~{sep='|' r2_files}".split("|")

ip_r1, ip_r2, ctl_r1, ctl_r2 = [], [], [], []
ip_pe, ctl_pe = None, None

for i, s in enumerate(samples):
    r1 = r1_paths[i]
    r2 = r2_paths[i]
    if s["is_control"]:
        ctl_r1.append(r1)
        ctl_r2.append(r2)
        ctl_pe = s["paired_end"]
    else:
        ip_r1.append(r1)
        ip_r2.append(r2)
        ip_pe = s["paired_end"]

def write_lines(fname, lines):
    with open(fname, "w") as f:
        f.write("\n".join(str(l) for l in lines))

write_lines("ip_r1.txt",  ip_r1)
write_lines("ip_r2.txt",  ip_r2)
write_lines("ctl_r1.txt", ctl_r1)
write_lines("ctl_r2.txt", ctl_r2)

with open("ip_pe.txt",  "w") as f: f.write(str(ip_pe).lower()  if ip_pe  is not None else "false")
with open("ctl_pe.txt", "w") as f: f.write(str(ctl_pe).lower() if ctl_pe is not None else "false")
PYEOF
    >>>

    output {
        Array[File] ip_r1          = read_lines("ip_r1.txt")
        Array[File] ip_r2          = read_lines("ip_r2.txt")
        Array[File] ctl_r1         = read_lines("ctl_r1.txt")
        Array[File] ctl_r2         = read_lines("ctl_r2.txt")
        Boolean     ip_paired_end  = read_boolean("ip_pe.txt")
        Boolean     ctl_paired_end = read_boolean("ctl_pe.txt")
    }

    runtime {
        cpu:         1
        memory:      "2 GB"
        disks:       "local-disk 10 SSD"
        docker:      docker
        singularity: singularity
        conda:       conda
    }
}
