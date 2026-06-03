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
#   2. Run with Caper or Terra alongside chip.wdl (imported below).
#   3. Outputs land in the same place as a normal chip.wdl run.
#
# The wrapper only handles download + handoff. All pipeline parameters
# (aligner, genome_tsv, peak caller, etc.) are passed through unchanged.
# ─────────────────────────────────────────────────────────────────────────────

https://raw.githubusercontent.com/ENCODE-DCC/chip-seq-pipeline2/refs/heads/master/chip.wdl as chip

struct SampleInfo {
    String accession       # SRR accession, e.g. "SRR10069059"
    Boolean paired_end     # true = PE, false = SE
    Boolean is_control     # true = input/IgG control, false = IP replicate
    String label           # human-readable label, used only for bookkeeping
}

workflow chip_sra_wrapper {

    meta {
        description: "Download FASTQs from SRA then run ENCODE ChIP-seq pipeline."
    }

    input {
        # ── Sample list ──────────────────────────────────────────────────────
        # Define 1 control + 1–4 IP replicates.
        # Order within IP replicates determines rep1, rep2, … assignment.
        Array[SampleInfo] samples

        # ── Genome / pipeline pass-through ───────────────────────────────────
        # These are forwarded directly to chip.wdl. At minimum supply
        # genome_tsv (easiest) or the individual genome files.
        File?   genome_tsv
        String  pipeline_type        # "tf" or "histone"
        String  title       = "Untitled"
        String  description = "No description"

        # Optional overrides (see chip.wdl for full parameter list)
        String  aligner              = "bowtie2"
        Boolean align_only           = false
        Boolean true_rep_only        = false
        Float   pval_thresh          = 0.01
        Float   idr_thresh           = 0.05

        # ── Runtime ──────────────────────────────────────────────────────────
        String  docker      = "encodedcc/chip-seq-pipeline:v2.2.2"
        String  singularity = "https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/chip-seq-pipeline_v2.2.2.sif"
        String  conda       = "encd-chip"

        # Resources for the download task
        Int     download_cpu         = 6
        Int     download_disk_gb     = 120   # per sample; SRA + FASTQ + headroom
        Int     download_mem_gb      = 8
    }

    # ── Download all samples in parallel ────────────────────────────────────
    scatter (s in samples) {
        call sra_to_fastq {
            input:
                accession    = s.accession,
                paired_end   = s.paired_end,
                cpu          = download_cpu,
                mem_gb       = download_mem_gb,
                disk_gb      = download_disk_gb,
                docker       = docker,
                singularity  = singularity,
                conda        = conda
        }
    }

    # ── Split downloaded FASTQs into IP replicates vs control ───────────────
    #
    # We collect R1 and R2 files from non-control samples as IP replicates
    # and from control samples as ctl inputs.  The chip.wdl expects individual
    # arrays per replicate (fastqs_rep1_R1, fastqs_rep2_R1, …), so we index
    # into the collected arrays below.
    #
    # Note: WDL 1.0 doesn't support filtering arrays natively, so we use
    # helper tasks for the split. For simplicity, control samples must be
    # listed LAST in the `samples` array (or use is_control flag – see task).

    call split_samples {
        input:
            accessions   = samples,
            r1_files     = sra_to_fastq.fastq_R1,
            r2_files     = sra_to_fastq.fastq_R2,
            docker       = docker,
            singularity  = singularity,
            conda        = conda
    }

    # ── Hand off to the ENCODE chip pipeline ────────────────────────────────
    #
    # Each IP replicate's R1 FASTQ goes into fastqs_repN_R1 as a single-
    # element array (one FASTQ file = one "technical replicate to merge").
    # The pipeline will merge if you supply multiple files per rep.
    #
    # We support up to 4 IP replicates here. Extend by adding more
    # fastqs_repN_R1 lines if needed.

    call chip.chip {
        input:
            # metadata
            title               = title,
            description         = description,
            pipeline_type       = pipeline_type,
            genome_tsv          = genome_tsv,

            # endedness: all replicates share the same setting for now
            # (extend SampleInfo and split_samples if you need mixed PE/SE)
            paired_end          = split_samples.ip_paired_end,
            ctl_paired_end      = split_samples.ctl_paired_end,

            # IP replicates – wrap each file in a single-element array
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

            # Control replicate (only first control used; pool if >1)
            ctl_fastqs_rep1_R1  = if length(split_samples.ctl_r1) > 0
                                    then [split_samples.ctl_r1[0]] else [],
            ctl_fastqs_rep1_R2  = if length(split_samples.ctl_r2) > 0
                                    then [split_samples.ctl_r2[0]] else [],

            # Pipeline flags
            aligner             = aligner,
            align_only          = align_only,
            true_rep_only       = true_rep_only,
            pval_thresh         = pval_thresh,
            idr_thresh          = idr_thresh,

            # Runtime
            docker              = docker,
            singularity         = singularity,
            conda               = conda
    }

    output {
        File    qc_report           = chip.report
        File    qc_json             = chip.qc_json
        Boolean qc_json_ref_match   = chip.qc_json_ref_match

        # Pass through the downloaded FASTQs so they land in the workspace
        # and can be reused in future runs without re-downloading.
        Array[File]  fastqs_R1      = sra_to_fastq.fastq_R1
        Array[File?] fastqs_R2      = sra_to_fastq.fastq_R2
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

        # prefetch caches the .sra file locally first (more reliable than
        # streaming directly with fasterq-dump on large files)
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

        # Compress outputs – chip.wdl expects .fastq.gz
        gzip -f *.fastq

        # Rename to predictable names regardless of SRA naming conventions
        # Expected outputs: {accession}_R1.fastq.gz and optionally _R2.fastq.gz
        if ls ~{accession}_1.fastq.gz 2>/dev/null; then
            mv ~{accession}_1.fastq.gz ~{accession}_R1.fastq.gz
        elif ls ~{accession}.fastq.gz 2>/dev/null; then
            # Single-end: fasterq-dump emits just accession.fastq
            mv ~{accession}.fastq.gz ~{accession}_R1.fastq.gz
        fi

        if ls ~{accession}_2.fastq.gz 2>/dev/null; then
            mv ~{accession}_2.fastq.gz ~{accession}_R2.fastq.gz
        fi

        echo "Done."
        ls -lh *.fastq.gz
    >>>

    output {
        File  fastq_R1 = "~{accession}_R1.fastq.gz"
        File? fastq_R2 = if paired_end
                            then "~{accession}_R2.fastq.gz"
                            else None
    }

    runtime {
        cpu:     cpu
        memory:  "~{mem_gb} GB"
        disks:   "local-disk ~{disk_gb} SSD"
        # sra-tools must be in the image; see README for a compatible image
        docker:      "ncbi/sra-tools:3.0.0"
        singularity: singularity
        conda:       "sra-tools"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task: split the flat scatter outputs into IP vs control arrays
#
# WDL 1.0 cannot filter arrays with inline expressions, so we do it in a
# tiny Python snippet and return the split lists via stdout files.
# ─────────────────────────────────────────────────────────────────────────────

task split_samples {
    input {
        Array[SampleInfo] accessions
        Array[File]       r1_files
        Array[File?]      r2_files
        String            docker
        String            singularity
        String            conda
    }

    # Serialise the sample metadata so Python can read it
    File accessions_json = write_json(accessions)

    command <<<
        python3 <<'PYEOF'
        import json, sys

        with open("~{accessions_json}") as f:
            samples = json.load(f)

        r1_paths = "~{sep='|' r1_files}".split("|")
        r2_raw   = "~{sep='|' select_all(r2_files)}"
        r2_paths = r2_raw.split("|") if r2_raw.strip() else []

        ip_r1, ip_r2, ctl_r1, ctl_r2 = [], [], [], []
        ip_pe, ctl_pe = None, None
        r2_idx = 0

        for i, s in enumerate(samples):
            r1 = r1_paths[i]
            # R2 only present when paired_end=true
            r2 = None
            if s["paired_end"]:
                r2 = r2_paths[r2_idx] if r2_idx < len(r2_paths) else None
                r2_idx += 1

            if s["is_control"]:
                ctl_r1.append(r1)
                if r2:
                    ctl_r2.append(r2)
                ctl_pe = s["paired_end"]
            else:
                ip_r1.append(r1)
                if r2:
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
        Array[File]  ip_r1         = read_lines("ip_r1.txt")
        Array[File]  ip_r2         = read_lines("ip_r2.txt")
        Array[File]  ctl_r1        = read_lines("ctl_r1.txt")
        Array[File]  ctl_r2        = read_lines("ctl_r2.txt")
        Boolean      ip_paired_end  = read_boolean("ip_pe.txt")
        Boolean      ctl_paired_end = read_boolean("ctl_pe.txt")
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
