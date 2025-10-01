set -x -e

######Tools
minimap2="/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411data/minimap2-2.30_x64-linux/minimap2"
meryl="/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/data/meryl/build/bin/meryl"
hifiasm="/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/data/hifiasm/hifiasm"
merqury="/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/data/merqury/merqury.sh"
######Directorys
work="/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/results"
scratch="/sc/arion/scratch/arayan01/projects/nanopore/2025-08-29_assemble_SEM_RS411/results"
rs411_ont_data="/sc/arion/projects/oscarlr/nanopore_data/data/08_26_2025_RS4-11_Ligation/RS4-11_4/20250826_1757_P2S-03054-B_PBE56067_1b60b28b/fastq_pass/"
sem_ont_data="/sc/arion/projects/oscarlr/nanopore_data/data/08_26_2025_SEM_Ligation/SEM_2-1/50826_1756_P2S-03054-A_PBE56017_762e1a7b/fastq_pass/"
hifi_data="/sc/arion/projects/oscarlr/jackie/IGH_epigenome/results/2025-04-28_assembling_SEM_RS4-11/run_hifiasm"
######Job parameters
jobs="${work}/jobs"
errmsg="${work}/errors" 
######Genomes 
ref="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/hg38_filtered.fa"
minimap2_ref="/sc/arion/scratch/arayan01/projects/personalized_epigenome/data/hg38_rfgen/hg38_filtered.mmi"
######

function samplelist {
    output="${work}/sample_list.txt"
    
    # Make sure directory exists
    mkdir -p "$(dirname "$output")"
    
    # Clear file before writing
    > "$output"
    
    # RS4-11 sample
    rs411_ont_data="${rs411_ont_data%/}"  # remove trailing slash
    rs411_fastqs=("${rs411_ont_data}"/*.fastq.gz)
    echo -e "RS4-11\tont\t${rs411_fastqs[*]}" >> "$output"
    
    # SEM sample
    sem_ont_data="${sem_ont_data%/}"  # remove trailing slash
    sem_fastqs=("${sem_ont_data}"/*.fastq.gz)
    echo -e "SEM\tont\t${sem_fastqs[*]}" >> "$output"
    
    echo "Sample list written to $output"

}

function index_ref {
   bsub -P acc_oscarlr -q premium -n 4 -W 4:00 -R 'rusage[mem=16000]' \
        -o $jobs/minimap_index.log -eo $errmsg/minimap_index.log \
        bash -c "$minimap2 -d $minimap2_ref $ref"

}
function align_minimap2 {
    input="${work}/sample_list.txt"
    rm -rf "${scratch}/alignments"
    output_dir="${scratch}/alignments"
    mkdir -p "$output_dir"

    while read -r sample_id type fastq_files; do
        output_bam="${output_dir}/${sample_id}.bam"
        bsub -P acc_oscarlr -q premium -n 10 -W 72:00 -R "rusage[mem=40000] span[hosts=1]" \
            -o "$jobs/${sample_id}_minimap2.log" \
            -eo "$errmsg/${sample_id}_minimap2.err" \
            bash -c "$minimap2 -ax map-ont $minimap2_ref $fastq_files | samtools sort -o $output_bam && samtools index $output_bam"
    done < "$input"
}
function bam_stats_ont {
    bam_dir="/sc/arion/scratch/arayan01/projects/nanopore/results/alignments"
    out_dir="${bam_dir}/stats_coverage"
    mkdir -p "$out_dir"

    # Write regions to a BED file for bedcov
    bed_file="${out_dir}/IG_regions.bed"
    cat > "$bed_file" <<EOF
chr14    105526583    106879811    IGH
chr22    21957025     23125032     IGL
chr2     88691673     90370929     IGK
EOF

    for bam in "$bam_dir"/*.bam; do
        sample=$(basename "$bam" .bam)
        stats_out="${out_dir}/${sample}_stats.txt"
        cov_out="${out_dir}/${sample}_coverage.txt"

        bsub -P acc_oscarlr -q premium -n 10 -W 72:00 -R "rusage[mem=40000] span[hosts=1]" \
             -o "${jobs}/${sample}.log" \
             -eo "${errmsg}/${sample}.err" \
             bash -c "
                 samtools stats $bam > $stats_out 
             "
    done
}
function bam_coverage_ont {
    bam_dir="/sc/arion/scratch/arayan01/projects/nanopore/results/alignments"
    out_dir="${bam_dir}/coverage"
    mkdir -p "$out_dir"

    # BED file (0-based starts)
    bed_file="${out_dir}/IG_regions.bed"
    regions=("IGH" "IGL" "IGK")
    printf "chr14\t105526582\t106879811\nchr22\t21957024\t23125032\nchr2\t88691672\t90370929\n" > "$bed_file"

    for bam in "$bam_dir"/*.bam; do
        sample=$(basename "$bam" .bam)
        out_file="${out_dir}/${sample}_coverage.txt"

        bsub -P acc_oscarlr -q premium -n 10 -W 72:00 -R "rusage[mem=40000] span[hosts=1]" \
             -o "${jobs}/${sample}_cov.log" \
             -eo "${errmsg}/${sample}_cov.err" \
             bash -c '
                 regions=("IGH" "IGL" "IGK")
                 bedcov_out=$(samtools bedcov '"$bed_file"' '"$bam"')
                 echo -e "Region\tChrom\tStart\tEnd\tRegionLength\tTotalBases\tMeanCoverage" > '"$out_file"'
                 i=0
                 while read -r chrom start end total; do
                     region_name=${regions[i]}
                     length=$((end-start))
                     mean=$(awk -v t=$total -v l=$length '\''BEGIN{printf "%.2f", t/l}'\'')
                     echo -e "${region_name}\t${chrom}\t${start}\t${end}\t${length}\t${total}\t${mean}" >> '"$out_file"'
                     ((i++))
                 done <<< "$bedcov_out"
             '
    done
}
function assemble_hifiasm_ont {
    sample_file="${work}/sample_list.txt"
    out_dir="${scratch}/ont_hifiasm_assemblies"
    mkdir -p "$out_dir"

    while read -r sample type fastqs; do
        # Only assemble ONT samples
        if [[ "$type" != "ont" ]]; then
            continue
        fi

        # Output prefix
        asm_out="${out_dir}/${sample}.asm"

        # Submit LSF job using --ont mode for noisy ONT reads
        bsub -P acc_oscarlr -q premium -n 10 -W 72:00 -R "rusage[mem=40000] span[hosts=1]" \
            -o "${jobs}/${sample}_hifiasm.log" \
            -eo "${errmsg}/${sample}_hifiasm.err" \
            bash -c "
                ${hifiasm} --ont -o ${asm_out} -t 10 ${fastqs} &&
                gfatools gfa2fa ${asm_out}.bp.hap1.p_ctg.gfa > ${asm_out}-hap1.fasta &&
                gfatools gfa2fa ${asm_out}.bp.hap2.p_ctg.gfa > ${asm_out}-hap2.fasta
            "
    done < "$sample_file"
}

#########################################
#########################################
#PACBIO DATA#############################
#########################################
#########################################

hifi_files="/sc/arion/scratch/arayan01/projects/nanopore/2025-08-29_assemble_SEM_RS411/data/Pacbio/hifi-reads"
SEM_WGS_INPUT="${hifi_files}/SEM-HiFi-WGS.fastq.gz"
SEM_HiFi_INPUT="${hifi_files}/SEM-Targeted-HiFi.fastq.gz"
RS4_HiFi_INPUT="${hifi_files}/RS411-Targeted-HiFi.fastq.gz"

hifi_data=(
    "${SEM_WGS_INPUT}"
    "${SEM_HiFi_INPUT}"
    "${RS4_HiFi_INPUT}"
)

# Define names for hifi reads
declare -A names
names["${SEM_WGS_INPUT}"]="SEM-HiFi-WGS"
names["${SEM_HiFi_INPUT}"]="SEM-Targeted-HiFi"
names["${RS4_HiFi_INPUT}"]="RS411-Targeted-HiFi"

assemblies=(
    "SEM-HiFi-WGS"
    "SEM-Targeted-HiFi"
    "RS411-Targeted-HiFi"
    "SEM-combined"
)

# map each asm name to its input FASTQ(s)
declare -A input_for=(
  ["SEM-HiFi-WGS"]="${SEM_WGS_INPUT}"
  ["SEM-Targeted-HiFi"]="${SEM_HiFi_INPUT}"
  ["RS411-Targeted-HiFi"]="${RS4_HiFi_INPUT}"
)

function run_hifiasm_pacbio {
    mkdir -p "${scratch}/pacbio_hifiasm"
    hifiasm_dir="${scratch}/pacbio_hifiasm"

    # Loop through each assembly
    for asm in "${assemblies[@]}"; do
        input_file="${input_for[$asm]}"   # Input FASTQ(s) for this assembly
        output_name="${asm}"              # Output folder name
        outdir="${hifiasm_dir}/${output_name}"

        mkdir -p "$outdir"

        if [[ "$asm" == "SEM-combined" ]]; then
            inputs="${SEM_WGS_INPUT} ${SEM_HiFi_INPUT}"
        else
            inputs="$input_file"
        fi

        # Submit each assembly as a separate job to LSF
        bsub -P  acc_oscarlr -q premium -n 10 -W 72:00 -R "rusage[mem=40000] span[hosts=1]" \
            -o "${jobs}/${asm}_pacbio_hifiasm.log" \
            -eo "${errmsg}/${asm}_pacbio_hifiasm.err" \
            bash -c "
                ${hifiasm} -o ${outdir}/${asm}.asm -t 10 ${inputs} &&
                gfatools gfa2fa ${outdir}/${asm}.asm.bp.hap1.p_ctg.gfa > ${outdir}/${output_name}-hap1.fasta &&
                gfatools gfa2fa ${outdir}/${asm}.asm.bp.hap2.p_ctg.gfa > ${outdir}/${output_name}-hap2.fasta
            "
    done
}

# Function to build meryl k-mer database from PacBio reads
function run_meryl_pacbio {
    for asm in "${assemblies[@]}"; do
        if [[ "$asm" == "SEM-combined" ]]; then
            continue
        fi

        # Assume input_for[$asm] contains both hap1 and hap2 fastq paths
        for hap in 1 2; do
            input_reads="${input_for[${asm}_hap${hap}]}"  # e.g., input_for["SEM-hap1"]
            outdir="${scratch}/pacbio_meryl_results/${asm}_hap${hap}"
            mkdir -p "$outdir/logs"

            # Remove old k-mer DB if exists
            if [[ -d "${outdir}/reads.meryl" ]]; then
                echo "Removing old meryl DB for $asm hap${hap}"
                rm -rf "${outdir}/reads.meryl"
            fi

            bsub -P acc_oscarlr -q premium -n 10 -W 24:00 -R "rusage[mem=40000] span[hosts=1]" \
                -o "${outdir}/logs/${asm}_hap${hap}_meryl.log" \
                -eo "${outdir}/logs/${asm}_hap${hap}_meryl.err" \
                bash -c "
                    export MERQURY=/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/data/merqury
                    export PATH=\$MERQURY/build/bin:\$PATH

                    ${meryl} count k=21 output ${outdir}/reads.meryl ${input_reads}
                "
        done
    done
}

# Function to run Merqury on PacBio hap1/hap2 assemblies (diploid mode)
function run_merqury_pacbio {
    local meryl_base="$scratch/pacbio_meryl_results"   # where reads.meryl lives
    local asm_base="$scratch/pacbio_hifiasm"           # where hap1/hap2 FASTAs live
    local out_base="$scratch/pacbio_merc_results"     # separate folder for Merqury outputs

    mkdir -p "$out_base"

    for asm in "${assemblies[@]}"; do
        if [[ "$asm" == "SEM-combined" ]]; then
            continue
        fi

        hap1_fasta="${asm_base}/${asm}/${asm}-hap1.fasta"
        hap2_fasta="${asm_base}/${asm}/${asm}-hap2.fasta"
        read_meryl="${meryl_base}/${asm}_hap1/reads.meryl"   # use hap1 reads.meryl
        outdir="${out_base}/${asm}"

        # Remove old Merqury output if it exists
        if [[ -d "${outdir}" ]]; then
            echo "Removing old Merqury output for $asm"
            rm -rf "${outdir:?}"
        fi

        # Create output + logs folder
        mkdir -p "${outdir}/logs"

        bsub -P acc_oscarlr -q premium -n 10 -W 24:00 -R "rusage[mem=40000] span[hosts=1]" \
            -o "${outdir}/logs/${asm}_merqury.log" \
            -eo "${outdir}/logs/${asm}_merqury.err" \
            bash -c "
                export MERQURY=/sc/arion/work/arayan01/project/nanopore/2025-08-29_assemble_SEM_RS411/data/merqury
                export PATH=\$MERQURY/build/bin:\$PATH

                # Run Merqury in diploid mode: reads.meryl + hap1 + hap2 FASTAs
                \$MERQURY/merqury.sh ${read_meryl} ${hap1_fasta} ${hap2_fasta} ${outdir}
            "
    done
}






function summarize_merqury_results {
    merqury_dir="${scratch}/merqury_results"
    output_table="${scratch}/merqury_summary.tsv"

    # Header
    echo -e "Assembly\tHaplotype\tQV\tCompleteness" > "$output_table"

    # Loop through assemblies
    for asm_dir in "$merqury_dir"/*; do
        asm_name=$(basename "$asm_dir")
        
        # Skip if not a directory
        [[ -d "$asm_dir" ]] || continue
        
        merqury_subdir="${asm_dir}/${asm_name}_merqury"
        
        for hap in hap1 hap2; do
            qv_file="${merqury_subdir}/${hap}.qv"
            if [[ -f "$qv_file" ]]; then
                qv=$(grep "QV" "$qv_file" | awk '{print $2}')
                comp=$(grep "completeness" "$qv_file" | awk '{print $2}')
                echo -e "${asm_name}\t${hap}\t${qv}\t${comp}" >> "$output_table"
            fi
        done
    done

    echo "Merqury summary table written to $output_table"
}








#ONT DATA#############################
#########################################
#samplelist
#index_ref
#align_minimap2
#bam_stats_ont
#bam_coverage_ont
#assemble_hifiasm_ont
#########################################
#PACBIO DATA#############################
#########################################
#run_hifiasm_pacbio 
#run_meryl_pacbio
run_merqury_pacbio
#summarize_merqury_results