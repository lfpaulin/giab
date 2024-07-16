#!/bin/bash
#SBATCH --job-name=HG008
#SBATCH --ntasks=24
#SBATCH --mem=48G
#SBATCH --output=log_hg008_2.out
#SBATCH --error=log_hg008_2.err

# CONDA
conda activate smaht

# tmp dir
TMPDIR=$1
# final dir
RESDIR=$2
# reference genome dir
REFERENCES=$3

# DATA
hg008_t_ont_ul="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/UCSC_ONT-UL_20231207/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_20231031.bam"
hg008_t_ont_ul_name="HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado0.4.3_sup4.2.0_5mCG_5hmCG_54x_UCSC_20231031.bam"
hg008_t_ont_ul_link="HG008_T-GRCh38_GIABv3-ONT_UL_54x.bam"


hg008_t_ont_wg="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/UCSC_ONT_20231003/HG008-T_GRCh38_GIABv3_ONT-R10.4.1-doradov0.3.4-5mCG-5hmC-63x_UCSC_20230905.bam"
hg008_t_ont_wg_name="HG008-T_GRCh38_GIABv3_ONT-R10.4.1-doradov0.3.4-5mCG-5hmC-63x_UCSC_20230905.bam"
hg008_t_ont_wg_link="HG008_T-GRCh38_GIABv3-ONT_WG_63x.bam"

hg008_nd_ont_wg="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x.bam"
hg008_nd_ont_wg_name="HG008-N-D_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_94x.bam"
hg008_nd_ont_wg_link="HG008_ND-GRCh38_GIABv3-ONT_94x.bam"

hg008_np_ont_wg="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam"
hg008_np_ont_wg_name="HG008-N-P_GRCh38-GIABv3_ONT-R1041-dorado_0.5.3_5mC_5hmC_41x.bam"
hg008_np_ont_wg_link="HG008_NP-GRCh38_GIABv3-ONT_41x.bam"


# REF GRCh38 GIABv3
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
REF38_GIABv3="${REFERENCES}/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"


# download data
data_download=(
    "${hg008_t_ont_ul},${hg008_t_ont_ul_name},${hg008_t_ont_ul_link}"
    "${hg008_t_ont_wg},${hg008_t_ont_wg_name},${hg008_t_ont_wg_link}"
    "${hg008_nd_ont_wg},${hg008_nd_ont_wg_name},${hg008_nd_ont_wg_link}"
    "${hg008_np_ont_wg},${hg008_np_ont_wg_name},${hg008_np_ont_wg_link}"
)


cd ${TMPDIR}
for datadl in ${data_download[@]};
do
    bam="$(echo ${datadl} | cut -d "," -f 1)"
    bai="${bam}.bai"
    name="$(echo ${datadl} | cut -d "," -f 2)"
    link="$(echo ${datadl} | cut -d "," -f 3)"
    # bam
    echo ${bam};
    wget ${bam};
    # name
    echo "${bai}";
    wget "${bai}";
    # links
    ln -T ${name} -s ${link}
    ln -T ${name}.bai -s ${link}.bai
done


# MERGE DATA
# THREADS = 24
hg008_t_ont_merge="HG008_T_merged-GRCh38_GIABv3-ONT_117x.bam"
samtools merge --threads 22 -o ${hg008_t_ont_merge}  ${hg008_t_ont_ul_link}  ${hg008_t_ont_wg_link} 
samtools index --threads 23 ${hg008_t_ont_merge}

hg008_n_ont_merge="HG008_N_merged-GRCh38_GIABv3-ONT_135x.bam"
samtools merge --threads 22 -o ${hg008_n_ont_merge}  ${hg008_nd_ont_wg_link}  ${hg008_np_ont_wg_link}
samtools index --threads 23 ${hg008_n_ont_merge}


# Analysis data
MODEL_GUPPY6="~/data/guppy_models/clair3_model-guppy_6" # not provided

giab_chunks_path="chunks"
chunk_list=(
    "chunk_001_chr1_1_100000000,${giab_chunks_path}/chunk_001_chr1_1_100000000.bed"
    "chunk_002_chr1_100000001_200000000,${giab_chunks_path}/chunk_002_chr1_100000001_200000000.bed"
    "chunk_003_chr1_200000001_248956422,${giab_chunks_path}/chunk_003_chr1_200000001_248956422.bed"
    "chunk_004_chr2_1_100000000,${giab_chunks_path}/chunk_004_chr2_1_100000000.bed"
    "chunk_005_chr2_100000001_200000000,${giab_chunks_path}/chunk_005_chr2_100000001_200000000.bed"
    "chunk_006_chr2_200000001_242193529,${giab_chunks_path}/chunk_006_chr2_200000001_242193529.bed"
    "chunk_007_chr3_1_100000000,${giab_chunks_path}/chunk_007_chr3_1_100000000.bed"
    "chunk_008_chr3_100000001_198295559,${giab_chunks_path}/chunk_008_chr3_100000001_198295559.bed"
    "chunk_009_chr4_1_100000000,${giab_chunks_path}/chunk_009_chr4_1_100000000.bed"
    "chunk_010_chr4_100000001_190214555,${giab_chunks_path}/chunk_010_chr4_100000001_190214555.bed"
    "chunk_011_chr5_1_100000000,${giab_chunks_path}/chunk_011_chr5_1_100000000.bed"
    "chunk_012_chr5_100000001_181538259,${giab_chunks_path}/chunk_012_chr5_100000001_181538259.bed"
    "chunk_013_chr6_1_100000000,${giab_chunks_path}/chunk_013_chr6_1_100000000.bed"
    "chunk_014_chr6_100000001_170805979,${giab_chunks_path}/chunk_014_chr6_100000001_170805979.bed"
    "chunk_015_chr7_1_100000000,${giab_chunks_path}/chunk_015_chr7_1_100000000.bed"
    "chunk_016_chr7_100000001_159345973,${giab_chunks_path}/chunk_016_chr7_100000001_159345973.bed"
    "chunk_017_chr8_1_100000000,${giab_chunks_path}/chunk_017_chr8_1_100000000.bed"
    "chunk_018_chr8_100000001_145138636,${giab_chunks_path}/chunk_018_chr8_100000001_145138636.bed"
    "chunk_019_chr9_1_100000000,${giab_chunks_path}/chunk_019_chr9_1_100000000.bed"
    "chunk_020_chr9_100000001_138394717,${giab_chunks_path}/chunk_020_chr9_100000001_138394717.bed"
    "chunk_021_chr10_1_100000000,${giab_chunks_path}/chunk_021_chr10_1_100000000.bed"
    "chunk_022_chr10_100000001_133797422,${giab_chunks_path}/chunk_022_chr10_100000001_133797422.bed"
    "chunk_023_chr11_1_100000000,${giab_chunks_path}/chunk_023_chr11_1_100000000.bed"
    "chunk_024_chr11_100000001_135086622,${giab_chunks_path}/chunk_024_chr11_100000001_135086622.bed"
    "chunk_025_chr12_1_100000000,${giab_chunks_path}/chunk_025_chr12_1_100000000.bed"
    "chunk_026_chr12_100000001_133275309,${giab_chunks_path}/chunk_026_chr12_100000001_133275309.bed"
    "chunk_027_chr13_1_100000000,${giab_chunks_path}/chunk_027_chr13_1_100000000.bed"
    "chunk_028_chr13_100000001_114364328,${giab_chunks_path}/chunk_028_chr13_100000001_114364328.bed"
    "chunk_029_chr14_1_100000000,${giab_chunks_path}/chunk_029_chr14_1_100000000.bed"
    "chunk_030_chr14_100000001_107043718,${giab_chunks_path}/chunk_030_chr14_100000001_107043718.bed"
    "chunk_031_chr15_1_100000000,${giab_chunks_path}/chunk_031_chr15_1_100000000.bed"
    "chunk_032_chr15_100000001_101991189,${giab_chunks_path}/chunk_032_chr15_100000001_101991189.bed"
    "chunk_033_chr16_1_90338345,${giab_chunks_path}/chunk_033_chr16_1_90338345.bed"
    "chunk_034_chr17_1_83257441,${giab_chunks_path}/chunk_034_chr17_1_83257441.bed"
    "chunk_035_chr18_1_80373285,${giab_chunks_path}/chunk_035_chr18_1_80373285.bed"
    "chunk_036_chr19_1_58617616,${giab_chunks_path}/chunk_036_chr19_1_58617616.bed"
    "chunk_037_chr20_1_64444167,${giab_chunks_path}/chunk_037_chr20_1_64444167.bed"
    "chunk_038_chr21_1_46709983,${giab_chunks_path}/chunk_038_chr21_1_46709983.bed"
    "chunk_039_chr22_1_50818468,${giab_chunks_path}/chunk_039_chr22_1_50818468.bed"
    "chunk_040_chrX_1_100000000,${giab_chunks_path}/chunk_040_chrX_1_100000000.bed"
    "chunk_041_chrX_100000001_156040895,${giab_chunks_path}/chunk_041_chrX_100000001_156040895.bed"
    "chunk_042_chrY_1_57227415,${giab_chunks_path}/chunk_042_chrY_1_57227415.bed"
    "chunk_043_chr16_KI270728v1_random_1_1872759,${giab_chunks_path}/chunk_043_chr16_KI270728v1_random_1_1872759.bed"
)


# ANALYSIS
data_analysis=(
    "${TMPDIR}/${hg008_t_ont_ul_link}"
    "${TMPDIR}/${hg008_t_ont_wg_link}"
    "${TMPDIR}/${hg008_nd_ont_wg_link}"
    "${TMPDIR}/${hg008_np_ont_wg_link}"
    "${TMPDIR}/${hg008_t_ont_merge}"
    "${TMPDIR}/${hg008_n_ont_merge}"
)


THREADS=24
cd ${RESDIR}
for hg008BAM in ${data_analysis[@]};
do
    # Sample ID
    sid=$(basename ${hg008BAM} | cut -d "." -f 1)
    # #################################################### #
    # Coverage: mosdepth
    echo "coverage: mosdepth"
    cd ${RESDIR}
    COVDIR="${RESDIR}/coverage"
    if [[ ! -d "${COVDIR}" ]];
    then
        mkdir ${COVDIR}
    fi
    cd ${COVDIR}
    WINSIZE="1000"
    mosdepth \
        --by ${WINSIZE} \
        --threads ${THREADS} \
        --no-per-base \
        --mapq 10 \
        ${sid} \
        ${hg008BAM}
    # #################################################### #
    # SNV clair3
    echo "SNV: clair3"
    cd ${RESDIR}
    SNVDIR="${RESDIR}/snv"
    if [[ ! -d "${SNVDIR}" ]];
    then
        mkdir ${SNVDIR}
    fi
    cd ${SNVDIR}
    # CONDA clair3
    conda activate clair
    list_chunk_merge="list_chunks_${sid}.txt"
    for chunk in ${chunk_list[@]};
    do
        chunk_name=$(echo ${chunk} | cut -d "," -f 1)
        chunk_path=$(echo ${chunk} | cut -d "," -f 2)
        output_dir="snv_${sid}_${chunk_name}"
        echo "${SNVDIR}/${output_dir}/merge_output.vcf.gz" >> "${list_chunk_merge}"
        run_clair3.sh \
            --ref_fn="${REF38_GIABv3}" \
            --model_path="${MODEL_GUPPY6}" \
            --threads="${THREADS}" \
            --platform="ont" \
            --bam_fn="${hg008BAM}" \
            --output="${output_dir}" \
            --bed_fn="${chunk_path}"
    done
    # SNV merge
    conda activate svwork
    SNV_OUTPUT="${sid}_snv.vcf.gz"
    # THREADS = 32
    bcftools concat \
        --allow-overlaps \
        --remove-duplicates \
        --file-list "${list_chunk_merge}" \
        --threads ${THREADS} | \
          bcftools view --apply-filters PASS --type snps --include "AF > 0.1" | \
          bcftools sort --max-mem 16G --output  ${SNVDIR}/${SNV_OUTPUT} --output-type z

    bcftools index --tbi --threads --threads ${THREADS} ${SNVDIR}/${SNV_OUTPUT}
    # #################################################### #
    SV sniffles
    echo "SV: sniffles"
    cd ${RESDIR}
    SVDIR="${RESDIR}/sv"
    if [[ ! -d "${SVDIR}" ]];
    then
        mkdir ${SVDIR}
    fi
    cd ${SVDIR}
    sniffles \
        --input "${hg008BAM}" \
        --vcf "${sid}.vcf.gz"  \
        --snf "${sid}.snf"  \
        --reference "${REF38_GIABv3}" \
        --output-rnames \
        --threads ${THREADS} \
        --sample-id "${sid}"
    # #################################################### #
done

# CNV data has pre-requisited from previous analysis
# CNV/LOH  spectre
# included in the spectre repo
REF_META=$4
BLACKLIST=$5

echo "CNV: spectre"
cd ${RESDIR}
CNVDIR="${RESDIR}/cnv"
if [[ ! -d "${CNVDIR}" ]];
then
    mkdir ${CNVDIR}
fi
cd ${CNVDIR}
data_analysis=(
    "${TMPDIR}/${hg008_t_ont_ul_link},${TMPDIR}/${hg008_n_ont_merge}"
    "${TMPDIR}/${hg008_t_ont_wg_link},${TMPDIR}/${hg008_n_ont_merge}"
    "${TMPDIR}/${hg008_t_ont_merge},${TMPDIR}/${hg008_n_ont_merge}"
)
for TN_SAMPLES in ${data_analysis[@]};
do
    # Control
    CONTROL_S=$(echo ${TN_SAMPLES} | cut -d "," -f 2)
    CONTROL_SAMPLE_ID=$(basename ${CONTROL_S} | cut -d "." -f 1)
    CONTROL_COV_FILE="../coverage/${CONTROL_SAMPLE_ID}.regions.bed.gz"
    CONTROL_SNV_FILE="../snv/${CONTROL_SAMPLE_ID}.vcf.gz"
    OUTPUT_DIR="cnv_${CONTROL_SAMPLE_ID}"
    if [[ ! -d "${OUTPUT_DIR}" ]];
    then
        mkdir ${OUTPUT_DIR}
    fi
    echo "${CONTROL_SAMPLE_ID}"
    spectre CNVCaller \
        --sample-id    "${CONTROL_SAMPLE_ID}"  \
        --coverage     "${CONTROL_COV_FILE}"   \
        --output-dir   "${OUTPUT_DIR}" \
        --reference    "${REFERENCES}"  \
        --metadata     "${REF_META}"   \
        --blacklist    "${BLACKLIST}"  \
        --snv          "${CONTROL_SNV_FILE}"   \
        --dev    2> ${OUTPUT_DIR}/log_cnv_240624_${SID}.txt

    # Tumor
    TUMOR_S=$(echo ${TN_SAMPLES} | cut -d "," -f 1)
    TUMOR_SAMPLE_ID=$(basename ${TUMOR_S} | cut -d "." -f 1)
    TUMOR_COV_FILE="../coverage/${TUMOR_SAMPLE_ID}.regions.bed.gz"
    TUMOR_SNV_FILE="../snv/${TUMOR_SAMPLE_ID}.vcf.gz"
    OUTPUT_DIR="cnv_${TUMOR_SAMPLE_ID}"
    CANCER="${CONTROL_SNV_FILE},${CONTROL_COV_FILE}"
    OUTPUT_DIR="cnv_${SID}"
    if [[ ! -d "${OUTPUT_DIR}" ]];
    then
        mkdir ${OUTPUT_DIR}
    fi
    echo "${SAMPLE_ID}"
    spectre CNVCaller \
        --sample-id    "${TUMOR_SAMPLE_ID}"  \
        --coverage     "${TUMOR_COV_FILE}"   \
        --output-dir   "${OUTPUT_DIR}" \
        --reference    "${FASTA_REF}"  \
        --metadata     "${REF_META}"   \
        --blacklist    "${BLACKLIST}"  \
        --snv          "${TUMOR_SNV_FILE}"   \
        --cancer       "${CANCER}"     \
        --dev    2> ${OUTPUT_DIR}/log_cnv_240626_${SID}.txt
done
