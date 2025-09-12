set -x -e

######Tools
minimap2="/sc/arion/work/arayan01/project/nanopore/data/minimap2-2.30_x64-linux/minimap2"
meryl="/sc/arion/work/arayan01/project/nanopore/data/meryl/build/bin/meryl"
hifiasm="/sc/arion/work/arayan01/project/nanopore/data/hifiasm/hifiasm"
######Directorys
work="/sc/arion/work/arayan01/project/nanopore/results"
scratch="/sc/arion/scratch/arayan01/projects/nanopore/results"
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




#samplelist
#index_ref
align_minimap2


#    # HiFi: SEM_RS4-11
#    for f in "${hifi_data}"/*.fastq.gz; do
#        [ -e "$f" ] || continue
#        echo -e "SEM_RS4-11\thifi\t$f" >> "$output"
#    done
    