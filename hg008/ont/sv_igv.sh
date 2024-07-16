#!/bin/bash
# tmp dir
TMPDIR=$1

# final dir
RESDIR=$2

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
REF38_GIABv3="/stornext/snfs4/next-gen/scratch/luis/references/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz"


# download data
data_download=(
    "${hg008_t_ont_ul},${hg008_t_ont_ul_name},${hg008_t_ont_ul_link}"
    "${hg008_t_ont_wg},${hg008_t_ont_wg_name},${hg008_t_ont_wg_link}"
    "${hg008_nd_ont_wg},${hg008_nd_ont_wg_name},${hg008_nd_ont_wg_link}"
    "${hg008_np_ont_wg},${hg008_np_ont_wg_name},${hg008_np_ont_wg_link}"
)


cd ${TMPDIR}
BED_USE="/stornext/snfs4/next-gen/scratch/luis/smaht/analysis/giab/hg008/240522/igv/extrat_reads_hg008.bed"
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
    
    BAM_OUT="$(echo ${link} | cut -d "." -f 1)_IGV"

    # UPDATE THE NAME OF THE FILE ...
    BAM_IN=${link}

    # add header
    samtools view -H ${BAM_IN} > tmp.sam;
    # add every single candidate from the bed
    cat ${BED_USE} | perl -lane 'print "$F[0]:$F[1]-$F[2]"' | xargs -n1 -t -I{} samtools view ${BAM_IN} {} >> tmp.sam;

    # make it bam and sort it
    samtools view -bSh tmp.sam | samtools sort -o ${BAM_OUT} -
    samtools index ${BAM_OUT}

    mv ${BAM_OUT}     ${RESDIR}/${BAM_OUT}
    mv ${BAM_OUT}.bai ${RESDIR}/${BAM_OUT}.bai

done
