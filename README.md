# *Brugia malayi*  Meta-Analysis
An in-depth analysis of available *B. malayi* RNA-Sequencing data. The files groups.txt and end_type_and_strandness.txt for each project are available in input_files

[![DOI](https://zenodo.org/badge/745300546.svg)](https://zenodo.org/doi/10.5281/zenodo.13694246)

## Table of Contents
* [List of Projects and Papers](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#list-of-projects-and-papers)
* [Download Reads](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#download-reads)
* [Preprocess Reads](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#preprocess-reads)
* [Confirm Strandness Using Salmon](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#confirm-strandness-using-salmon)
* [Align Reads to Reference Genome](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#align-preprocessed-reads-to-reference-genome)
* [Merge Bam Files](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#merge-bam-file----required-for-chung_lifecycle_prjna294263)
* [Generate Counts](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#generate-counts)
* [R Session Info](https://github.com/christopher-holt/bmalayi_meta_analysis?tab=readme-ov-file#r-session-info)

## List of Projects and Papers
| Project ID | Paper | SRA BioProject |
| --------| ----- | -------------- |
| choi_PRJEB2709| A Deep Sequencing Approach to Comparatively Analyze the Transcriptome of Lifecycle Stages of the Filarial Worm, *Brugia malayi* | PRJEB2709 |
| ballesteros_PRJNA303987 | The Effects of Ivermectin on *Brugia malayi* Females In Vitro: A Transcriptomic Approach (Study I)| PRJNA303987 |
| ballesteros_PRJNA303986 | The Effects of Ivermectin on *Brugia malayi* Females In Vitro: A Transcriptomic Approach (Study II)| PRJNA303986 |
| ballesteros_PRJNA294426| The Effect of In Vitro Cultivation on the Transcriptome of Adult *Brugia malayi*| PRJNA294426 |
| libro_PRJNA329497| Characterization of innate immunity genes in the parasitic nematode *Brugia malayi*| PRJNA329497 |
| grote_PRJNA344486| Defining *Brugia malayi* and Wolbachia symbiosis by stage-specific dual RNA-Seq| PRJNA344486 |
| maclean_PRJNA388112| Effects of diethylcarbamazine and ivermectin treatment on *Brugia malayi* gene expression in infected gerbils (*Meriones unguiculatus*)| PRJNA388112 |
| chung_PRJNA294263| Drug Repurposing of Bromodomain Inhibitors as Potential Novel Therapeutic Leads for Lymphatic Filariasis Guided by Multispecies Transcriptomics| PRJNA294263 |
| grote_PRJNA557263| Prediction pipeline for discovery of regulatory motifs associated with *Brugia malayi* molting | PRJNA557263 |
| chevignon_PRJEB40568| Dual RNAseq analyses at soma and germline levels reveal evolutionary innovations in the elephantiasis-agent *Brugia malayi*, and adaptation of its *Wolbachia* endosymbionts| PRJEB40568 |
| quek_PRJNA772674 | *Wolbachia* depletion blocks transmission of lymphatic filariasis by preventing chitinase-dependent parasite exsheathment| PRJNA772674 |
| airs_PRJNA548881 | Spatial transcriptomics reveals antiparasitic targets associated with essential behaviors in the human parasite *Brugia malayi* | PRJNA548881 |




## Download Reads
```bash
PROJECT_ID=
WORKING_DIR=
PACKAGE_DIR=

SRA_TOOLKIT_BIN_DIR=${PACKAGE_DIR}/sratoolkit-2.10.9/bin/
READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}
SRR_LIST=${READS_DIR}/groups.txt

THREADS=1

mkdir -p ${OUTPUT_DIR}

for SRR in $(cat ${SRR_LIST}); do
    SRR_ID=$(echo ${SRR} | cut -d',' -f1)

    echo -e "
    ${SRA_TOOLKIT_BIN_DIR}/prefetch ${SRR_ID} --max-size 100G --output-directory ${READS_DIR}
    " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${READS_DIR} -N prefetch.${SRR_ID}

    echo -e "
    ${SRA_TOOLKIT_BIN_DIR}/fastq-dump --split-files ${SRR_ID} --gzip -O ${READS_DIR}/${SRR_ID}
    " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${READS_DIR} -N fastq.dump.${SRR_ID} -hold_jid prefetch.${SRR_ID}

done

```

## Preprocess Reads

```bash
WORKING_DIR=
PROJECT_ID=
PACKAGE_DIR=

READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}
SRR_LIST=${READS_DIR}/groups.txt
FASTP_BIN_DIR=${PACKAGE_DIR}/fastp-0.22.0/bin

THREADS=4

mkdir -p ${READS_DIR}
end_type=$(cat ${READS_DIR}/end_type_and_strandness.txt | cut -d',' -f1)


if [ ${end_type} == "single" ]; then
    echo "Single End"
    ## Single End Reads
    for SRR in $(cat ${SRR_LIST}); do
        SRR_ID=$(echo ${SRR} | cut -d',' -f1)

        echo -e "
        ${FASTP_BIN_DIR}/fastp.0.22.0 -w ${THREADS} -i ${READS_DIR}/${SRR_ID}/${SRR_ID}_1.fastq.gz -o ${READS_DIR}/${SRR_ID}/${SRR_ID}_fastp_1.fastq.gz
        " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${READS_DIR} -N preprocess.${SRR_ID} -hold_jid fastq.dump.${SRR_ID}

    done

else

    ## Paired End Reads
    for SRR in $(cat ${SRR_LIST}); do
        SRR_ID=$(echo ${SRR} | cut -d',' -f1)


        echo -e "
        ${FASTP_BIN_DIR}/fastp.0.22.0 -w ${THREADS} -i ${READS_DIR}/${SRR_ID}/${SRR_ID}_1.fastq.gz -I ${READS_DIR}/${SRR_ID}/${SRR_ID}_2.fastq.gz -o ${READS_DIR}/${SRR_ID}/${SRR_ID}_fastp_1.fastq.gz -O ${READS_DIR}/${SRR_ID}/${SRR_ID}_fastp_2.fastq.gz
        " | qsub -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=30G -wd ${READS_DIR} -N preprocess.${SRR_ID} -hold_jid fastq.dump.${SRR_ID}

    done
fi
```

## Confirm Strandness Using Salmon
```bash
WORKING_DIR=
PROJECT_ID=
PACKAGE_DIR=
REFERENCE_DIR=


READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}
SRR_LIST=${READS_DIR}/groups.txt

SALMON_BIN_DIR=${PACKAGE_DIR}/salmon-0.13.1/bin/
SALMON_OUTPUT_DIR=${WORKING_DIR}/salmon_strandness/${PROJECT}
BMALAYI_TRANSCRIPTS=${REFERENCE_DIR}/brugia_malayi.PRJNA10729.WBPS18.mRNA_transcripts.fa
BMALAYI_TRANSCRIPTS_INDEX=${REFERENCE_DIR}/brugia_malayi/brugia_malayi.PRJNA10729.WBPS18.mRNA_transcripts_salmon_index
end_type=$(cat ${READS_DIR}/end_type_and_strandness.txt | cut -d',' -f1)

THREADS=4

## Index Transcriptome - Only have to run once
${SALMON_BIN_DIR}/salmon index -t ${BMALAYI_TRANSCRIPTS} -i ${BMALAYI_TRANSCRIPTS_INDEX}

mkdir -p ${SALMON_OUTPUT_DIR}

for SRR in $(cat ${SRR_LIST}); do
    SRR_ID=$(echo $SRR | cut -d',' -f1)
    SAMPLE_NAME=$(echo $SRR | cut -d',' -f2)

    INPUT_DIR=${READS_DIR}/${SRR_ID}

    ## Single End Not Stranded
    if [ ${end_type} == "single" ]; then

        FASTQ=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_1.fastq.gz)

        echo -e "
        ${SALMON_BIN_DIR}/salmon quant -i ${BMALAYI_TRANSCRIPTS_INDEX} -l A -r ${FASTQ} --validateMappings -o ${SALMON_OUTPUT_DIR}
        " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${SALMON_OUTPUT_DIR} -N salmon.${SAMPLE_NAME} -hold_jid preprocess.${SRR_ID}

        ## Paired end
    elif [ ${end_type} == "paired" ]; then
            FASTQ1=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_1.fastq.gz)
            FASTQ2=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_2.fastq.gz)

             echo -e "
             ${SALMON_BIN_DIR}/salmon quant -i ${BMALAYI_TRANSCRIPTS_INDEX} -l A -1 ${FASTQ1} -2 ${FASTQ2} --validateMappings -o ${SALMON_OUTPUT_DIR}
             " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${SALMON_OUTPUT_DIR} -N salmon.${SAMPLE_NAME} -hold_jid preprocess.${SRR_ID}
         fi
done
```


## Align Preprocessed Reads to Reference Genome
```bash
WORKING_DIR=
PROJECT_ID=
PACKAGE_DIR=
REFERENCE_DIR=

READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}/
SRR_LIST=${READS_DIR}/groups.txt
HISAT2_DIR=${PACKAGE_DIR}/hisat2-2.1.0
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
BAM_DIR=${WORKING_DIR}/bam/${PROJECT}
REFERENCE_FASTA=${REFERENCE_DIR}/brugia_malayi/b_malayi.AF538716.AE017321.PRJNA10729.WS276.genomic.fa
REFERENCE_FASTA_BASE=${REFERENCE_DIR}/indexed_reference/b_malayi.AF538716.AE017321.PRJNA10729.WS276.genomic
end_type=$(cat ${READS_DIR}/end_type_and_strandness.txt | cut -d',' -f1)
strandness=$(cat ${READS_DIR}/end_type_and_strandness.txt | cut -d',' -f2)

THREADS=4

## Create Index for Hisat2

${HISAT2_DIR}/hisat2-build ${REFERENCE_GENOME} ${REFERENCE_FASTA_BASE}

sort_and_index_bam_file(){
       echo -e " ${SAMTOOLS_BIN_DIR}/samtools sort -@ ${THREADS} -o ${BAM} && ${SAMTOOLS_BIN_DIR}/samtools index -@ ${THREADS} ${BAM}"
    }

mkdir -p ${BAM_DIR}

for SRR in $(cat ${SRR_LIST}); do
    SRR_ID=$(echo $SRR | cut -d',' -f1)
    SAMPLE_NAME=$(echo $SRR | cut -d',' -f2)

    INPUT_DIR=${READS_DIR}/${SRR_ID}

    if [[ ${PROJECT} = "chung_lifecycle_PRJNA294263" ]]; then
        BAM=${BAM_DIR}/${SRR_ID}.sorted.bam
    else
        BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam
    fi

    ## Single End Not Stranded
    if [ ${end_type} == "single" ]; then

        echo "Single End, Not Stranded"
        FASTQ=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_1.fastq.gz)

        echo -e "
        ${HISAT2_DIR}/hisat2 -p ${THREADS} --max-intronlen 50000 -x ${REFERENCE_FASTA_BASE} -U ${FASTQ} | $(sort_and_index_bam_file)
        " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${BAM_DIR} -N hisat2.align.${SAMPLE_NAME} -hold_jid preprocess.${SRR_ID}

        ## Paired end Reversely Stranded
    elif [ ${end_type} == "paired" ]; then
        FASTQ1=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_1.fastq.gz)
        FASTQ2=$(echo ${INPUT_DIR}/${SRR_ID}_fastp_2.fastq.gz)
        if [ ${strandness} == "RF" ]; then

            echo "Paired End, reversely stranded"


             echo -e "
             ${HISAT2_DIR}/hisat2 -p ${THREADS} --rna-strandness RF --max-intronlen 50000 -x ${REFERENCE_FASTA_BASE} -1 ${FASTQ1} -2 ${FASTQ2} | $(sort_and_index_bam_file)
             " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${BAM_DIR} -N hisat2.align.${SAMPLE_NAME} -hold_jid preprocess.${SRR_ID}



         elif [ ${strandness} == "none" ]; then

             echo "Paired End, Not Stranded"

             echo -e "
             ${HISAT2_DIR}/hisat2 -p ${THREADS} --max-intronlen 50000 -x "$REFERENCE_FASTA_BASE" -1 ${FASTQ1} -2 ${FASTQ2} | $(sort_and_index_bam_file)
             " | qsub -V -P jhotopp-gcid-proj4b-filariasis -pe thread ${THREADS} -q threaded.q -l mem_free=50G -wd ${BAM_DIR} -N hisat2.align.${SAMPLE_NAME} -hold_jid preprocess.${SRR_ID}
         fi

    fi


done
```
## Merge Bam File -- Required for chung_lifecycle_PRJNA294263
```bash
PROJECT_ID=chung_lifecycle_PRJNA294263
WORKING_DIR=
PACKAGE_DIR=

READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}/
SRR_LIST=${READS_DIR}/srr.id.list
SAMTOOLS_BIN_DIR=${PACKAGE_DIR}/samtools-1.9/bin
BAM_DIR=${WORKING_DIR}/${PROJECT}

THREADS=4



mkdir -p ${BAM_DIR}

for SAMPLE_NAME in $(cat ${SRR_LIST} | cut -d',' -f2 | sort | uniq ); do

    SRR_ID=$(cat ${SRR_LIST} | grep ${SAMPLE_NAME} | cut -d',' -f1)


    while IFS=' ' read -ra ADDR; do
        for i in "${ADDR[@]}"; do
            echo ${BAM_DIR}/${i}.sorted.bam >> ${BAM_DIR}/${SAMPLE_NAME}.txt
        done
    done <<< "$SRR_ID"

    BAM=${BAM_DIR}/${SAMPLE_NAME}.bam
    SORTED_BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam

	echo -e "
    "$SAMTOOLS_BIN_DIR"/samtools merge -f -@ ${THREADS} -b ${BAM_DIR}/${SAMPLE_NAME}.txt ${BAM}; \
        "$SAMTOOLS_BIN_DIR"/samtools sort -o ${SORTED_BAM} ${BAM}; \
        "$SAMTOOLS_BIN_DIR"/samtools index ${SORTED_BAM}; \
        " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -l mem_free=20G -N samtools_merge.${SAMPLE_NAME} -wd ${BAM_DIR}

done
```
## Generate Counts
```bash
WORKING_DIR=
PROJECT_ID=
PACKAGE_DIR=
REFERENCE_DIR=

PYTHON_BIN_DIR=${PACKAGE_DIR}/python-3.8.2/bin
READS_DIR=${WORKING_DIR}/reads/${PROJECT_ID}
SRR_LIST=${READS_DIR}/groups.txt
GFF3=${REFERENCE_DIR}/brugia_malayi/b_malayi.AF538716.AE017321.PRJNA10729.WS276.annotations.wormbase.gff3
BAM_DIR=${WORKING_DIR}/bam/${PROJECT_ID}
COUNTS_DIR=${WORKING_DIR}/counts/${PROJECT_ID}
OUTFILE=${COUNTS_DIR}/combined.final.counts
SRA_ID=$(cat ${SRR_LIST} | awk -F',' '{ print $NF}' | sort | uniq)
strandness=$(cat ${READS_DIR}/end_type_and_strandness.txt | cut -d',' -f2)

mkdir -p ${COUNTS_DIR}

THREADS=4

if [ ${strandness} == "RF" ]; then
    STRANDNESS=reverse
else
    STRANDNESS=no
fi

for SAMPLE_NAME in $(cat ${SRR_LIST} | cut -d',' -f2 | sort | uniq); do


    BAM=${BAM_DIR}/${SAMPLE_NAME}.sorted.bam

    echo -e "
    ${PYTHON_BIN_DIR}/htseq-count -n ${THREADS} -s ${STRANDNESS} --max-reads-in-buffer 3000000000 -r pos --nonunique none -f bam -m union -t gene --idattr ID ${BAM} ${GFF3} | awk -v a=${SAMPLE_NAME} '{print \$0, a}' | awk -v b=${SRA_ID} '{print \$0, b}' | sed -e 's/ /\t/g' >> ${OUTFILE}
    " | qsub -V -q threaded.q -pe thread ${THREADS} -P jhotopp-gcid-proj4b-filariasis -N htseq_counts.${SAMPLE_NAME} -wd ${COUNTS_DIR} -l mem_free=70G -hold_jid hisat2.align.${SAMPLE_NAME}

done
```

## R Session Info
```r
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 14.5

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] devtools_2.4.5           usethis_2.1.6            WGCNA_1.71              
 [4] fastcluster_1.2.3        dynamicTreeCut_1.63-1    variancePartition_1.26.0
 [7] BiocParallel_1.30.4      pvclust_2.2-0            forcats_0.5.2           
[10] stringr_1.4.1            dplyr_1.0.10             purrr_0.3.5             
[13] readr_2.1.3              tidyr_1.2.1              tibble_3.1.8            
[16] ggplot2_3.4.0            tidyverse_1.3.2          ggdendro_0.1.23         
[19] FactoMineR_2.6           edgeR_3.38.4             limma_3.52.4            

loaded via a namespace (and not attached):
  [1] utf8_1.2.2             tidyselect_1.2.0       lme4_1.1-31           
  [4] RSQLite_2.2.18         AnnotationDbi_1.58.0   htmlwidgets_1.5.4     
  [7] grid_4.2.1             munsell_0.5.0          preprocessCore_1.58.0 
 [10] codetools_0.2-18       interp_1.1-3           DT_0.26               
 [13] miniUI_0.1.1.1         withr_2.5.0            colorspace_2.0-3      
 [16] Biobase_2.56.0         knitr_1.40             rstudioapi_0.14       
 [19] leaps_3.1              stats4_4.2.1           Rdpack_2.4            
 [22] emmeans_1.8.2          GenomeInfoDbData_1.2.8 bit64_4.0.5           
 [25] coda_0.19-4.1          vctrs_0.5.1            generics_0.1.3        
 [28] xfun_0.34              timechange_0.1.1       R6_2.5.1              
 [31] doParallel_1.0.17      GenomeInfoDb_1.32.4    locfit_1.5-9.6        
 [34] bitops_1.0-7           cachem_1.0.6           assertthat_0.2.1      
 [37] vroom_1.6.0            promises_1.2.0.1       scales_1.2.1          
 [40] nnet_7.3-18            googlesheets4_1.0.1    gtable_0.3.1          
 [43] multcompView_0.1-8     processx_3.8.0         rlang_1.0.6           
 [46] scatterplot3d_0.3-42   GlobalOptions_0.1.2    splines_4.2.1         
 [49] impute_1.70.0          gargle_1.2.1           broom_1.0.1           
 [52] checkmate_2.1.0        yaml_2.3.6             reshape2_1.4.4        
 [55] modelr_0.1.10          backports_1.4.1        httpuv_1.6.6          
 [58] Hmisc_4.7-1            tools_4.2.1            ellipsis_0.3.2        
 [61] gplots_3.1.3           RColorBrewer_1.1-3     BiocGenerics_0.42.0   
 [64] sessioninfo_1.2.2      Rcpp_1.0.9             plyr_1.8.8            
 [67] base64enc_0.1-3        progress_1.2.2         zlibbioc_1.42.0       
 [70] RCurl_1.98-1.9         ps_1.7.2               prettyunits_1.1.1     
 [73] rpart_4.1.19           deldir_1.0-6           urlchecker_1.0.1      
 [76] S4Vectors_0.34.0       haven_2.5.1            ggrepel_0.9.2         
 [79] cluster_2.1.4          fs_1.5.2               magrittr_2.0.3        
 [82] data.table_1.14.4      circlize_0.4.15        reprex_2.0.2          
 [85] googledrive_2.0.0      mvtnorm_1.1-3          matrixStats_0.62.0    
 [88] pkgload_1.3.2          hms_1.1.2              mime_0.12             
 [91] evaluate_0.18          xtable_1.8-4           pbkrtest_0.5.1        
 [94] RhpcBLASctl_0.21-247.1 jpeg_0.1-9             readxl_1.4.1          
 [97] IRanges_2.30.1         gridExtra_2.3          shape_1.4.6           
[100] compiler_4.2.1         KernSmooth_2.23-20     crayon_1.5.2          
[103] minqa_1.2.5            htmltools_0.5.3        later_1.3.0           
[106] tzdb_0.3.0             Formula_1.2-4          lubridate_1.9.0       
[109] DBI_1.1.3              dbplyr_2.2.1           MASS_7.3-58.1         
[112] boot_1.3-28            Matrix_1.5-3           cli_3.4.1             
[115] rbibutils_2.2.10       parallel_4.2.1         pkgconfig_2.0.3       
[118] flashClust_1.01-2      foreign_0.8-83         xml2_1.3.3            
[121] foreach_1.5.2          XVector_0.36.0         estimability_1.4.1    
[124] rvest_1.0.3            callr_3.7.3            digest_0.6.30         
[127] Biostrings_2.64.1      rmarkdown_2.18         cellranger_1.1.0      
[130] htmlTable_2.4.1        curl_4.3.3             shiny_1.7.3           
[133] gtools_3.9.3           nloptr_2.0.3           lifecycle_1.0.3       
[136] nlme_3.1-160           jsonlite_1.8.3         aod_1.3.2             
[139] fansi_1.0.3            pillar_1.8.1           lattice_0.20-45       
[142] KEGGREST_1.36.3        fastmap_1.1.0          httr_1.4.4            
[145] pkgbuild_1.3.1         survival_3.4-0         GO.db_3.15.0          
[148] glue_1.6.2             remotes_2.4.2          png_0.1-7             
[151] iterators_1.0.14       bit_4.0.5              stringi_1.7.8         
[154] profvis_0.3.7          blob_1.2.3             latticeExtra_0.6-30   
[157] caTools_1.18.2         memoise_2.0.1    
```
