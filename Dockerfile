#################### SVE IMAGE #######################
FROM srp33/somatic_wgs_environment:0.1.0

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD SCRIPTS ########################
ADD src/system_commands/echo_usage.sh /usr/local/bin/echo_usage
ADD src/system_commands/usage_functions.sh /usr/local/bin/usage_functions
ADD src/system_commands/check_functions.sh /usr/local/bin/check_functions
ADD src/system_commands/entrypoint.sh /usr/local/bin/wgs/entrypoint
ADD src/system_commands/check_permissions.py /

################ ADD WGS SCRIPTS ######################
ADD src/wgs/bwa_mem_align.sh /usr/local/bin/wgs/bwa_mem_align
ADD src/wgs/call_somatic_variants_gatk4.sh /usr/local/bin/wgs/call_somatic_variants_gatk4
ADD src/wgs/call_somatic_variants_strelka.sh /usr/local/bin/wgs/call_somatic_variants_strelka
ADD src/wgs/index_bam.sh /usr/local/bin/wgs/index_bam
ADD src/wgs/mark_duplicates.sh /usr/local/bin/wgs/mark_duplicates
ADD src/wgs/merge_bams.sh /usr/local/bin/wgs/merge_bams
ADD src/wgs/slice_bam.sh /usr/local/bin/wgs/slice_bam
ADD src/wgs/sort_bam.sh /usr/local/bin/wgs/sort_bam
ADD src/wgs/add_read_groups.sh /usr/local/bin/wgs/add_read_groups
ADD src/wgs/samtools_mpileup.sh /usr/local/bin/wgs/samtools_mpileup
ADD src/wgs/base_recalibrator.sh /usr/local/bin/wgs/base_recalibrator
ADD src/wgs/apply_bqsr.sh /usr/local/bin/wgs/apply_bqsr
ADD src/wgs/call_structural_variants_manta.sh /usr/local/bin/wgs/call_structural_variants_manta
ADD src/wgs/samblast.sh /usr/local/bin/wgs/samblast
ADD src/wgs/call_structural_variants_lumpy.sh /usr/local/bin/wgs/call_structural_variants_lumpy
ADD src/wgs/call_structural_variants_delly.sh /usr/local/bin/wgs/call_structural_variants_delly
ADD src/wgs/svtype_vcf.sh /usr/local/bin/wgs/svtype_vcf
ADD src/wgs/run_survivor.sh /usr/local/bin/wgs/run_survivor

################## SETUP WORKDIR #######################
WORKDIR /data
RUN chmod 777 /data

######## OVERRIDE PARLIAMENT2 ENTRYPOINT ##############
#ENTRYPOINT ["entrypoint"]

##################### RUN BWA #########################
CMD echo_usage
