#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/hpcmeta
========================================================================================
 nf-core/hpcmeta Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/hpcmeta
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline with ONP reads is as follows:

    nextflow run -c custom.conf fullpath/Metagenome-Workflow/main.nf -profile uiuc_singularity -qs 3 -resume -with-report -with-trace


    Mandatory arguments:
      --reads                       Path to ONP data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
    Trimming Options:
      --skipTrim                    To enable the pipeline to run porechop trimming. Possible values: true, false. Default: true.
      
    Filtering Options:
      --minimapRemoveHost           To enable the pipeline to run this analysis with minimap. Possible values: true, false. Default: false.
      --host                        Path to host database to use for filtering reads. It should be in fasta format.
      --mergeSanizited              To enable vsearch merge. Possible values: true, false. Default: false.

    Kraken Options:      
      --runKraken2                  To enable the pipeline to run this analysis. Possible values: true, false. Default: false.
      --Kraken2DB                   It should be the path of a folder containing a kraken-formatted database

    Assembly Options:
      --runAssembly                 To enable the pipeline to run this analysis. Possible values: true, false. Default: false.
      --assembler                   Assembly tool. Possible values: metaflye.
      
    Diamond Options:
      --runDiamond                  To enable the pipeline to run this analysis. Possible values: true, false. Default: false. Must also enable assembly.
      --diamondDB                   Path to the Diamond database to use for taxonomic assignment. It should be diamond-indexed.
    
    Phylogenetic profiling and Taxonomic binning options
      --runMetaPhlan2               To enable the pipeline to run this analysis. Possible values: true, false. Default: false. Must also enable filtering.  
      --runMetaBAT                  To enable the pipeline to run this analysis. Possible values: true, false. Default: false. Must also enable assembly.   

    Blobtools Options:
      --runBlobtools                To enable the pipeline to run this analysis. Possible values: true, false. Default: false. Must also enable assembly,Diamond.
      --blobscript                  Path to daa_to_tagc.pl
      --blobDB                      Path to the Diamond database to use for taxonomic assignment. It should be diamond-indexed.
      --blobTaxIds                  Path to the tabular file with taxids that is companion to the blobDB.
      --blobphylum                  Path to the colors_phylum.txt file
      --blobsuperkingdom            Path to the colors_Superkingdom.txt file

    Other options:
      --singleEnd                   Specifies that the input is single end reads. Default false
      --inputType                   Specifies the type of input. Default "onp"
      --outdir                      Path of the output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail               Same as --email, except only send mail if the workflow is not successful
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
  
    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SANITY CHECKS FOR OPTIONAL STEPS
 */

if (params.runAssembly && params.assembler != 'metaflye') exit 1, "Unknown value for param assembler"
if (params.minimapRemoveHost && ! params.host) exit 1, "Must provide params host in combination with param minimapRemoveHost"
if (params.runKraken2 && ! params.Kraken2DB) exit 1, "Must provide params Kraken2DB in combination with param runKraken2"
if (params.runDiamond && ! params.diamondDB) exit 1, "Must provide params diamondDB in combination with param runDiamond"
if (! params.runDiamond && params.runBlobtools) exit 1, "Cannot skip Diamond in combination with run Blobtools"
if (params.runBlobtools && ! params.blobscript)  exit 1, "Must provide params blobscript in combination with run Blobtools"
if (params.runBlobtools && ! params.blobDB)      exit 1, "Must provide params blobDB in combination with run Blobtools"
if (params.runBlobtools && ! params.blobTaxIds)  exit 1, "Must provide params blobTaxIds in combination with run Blobtools"
if (params.runBlobtools && ! params.blobphylum)  exit 1, "Must provide params blobphylum in combination with run Blobtools"
if (params.runBlobtools && ! params.blobsuperkingdom)  exit 1, "Must provide params blobsuperkingdom in combination with run Blobtools"

// Check if genome exists in the config file
// if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
//     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
// }


// TODO nf-core: Add any reference files that are needed
// Configurable reference genomes
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the channel below in a process, define the following:
//   input:
//   file fasta from ch_fasta
//
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * Create a channel for input read files
 */

Channel
       .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
       .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
       .into { read_files_qc; read_files_trimming }


// take care of persistent files

params.genomesDir           = "${workflow.workDir}/HostGenome-nf-produced"
genomeStore                 = params.genomesDir

if (params.minimapRemoveHost &&  params.host) { 
   genomeFile                  = file(params.host)
   genomePrefix                = genomeFile.getBaseName()
}

if (params.runBlobtools) {
   blobDBFile                  = file(params.blobDB)
   blobTaxFile                 = file(params.blobTaxIds)
   blobPhylumFile              = file(params.blobphylum)
   blobKingdomFile             = file(params.blobsuperkingdom)
   blobDBPrefix                = blobDBFile.getBaseName()    
}
// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// nf-core: Report custom parameters here
summary['Reads']            = params.reads
summary['Input Type']       = params.inputType
summary['Remove Host']      = params.minimapRemoveHost
if (params.minimapRemoveHost &&  params.host) { summary['RefGenome Host']   = params.host }
summary['Run Kraken2']      = params.runKraken2
if (params.runKraken2 && params.Kraken2DB)    { summary['Kraken2 DB']   = params.Kraken2DB }
summary['Run Diamond']      = params.runDiamond
if (params.runDiamond && params.diamondDB)    { summary['Diamond DB']   = params.diamondDB }
summary['Run MetaPhaln2']   = params.runMetaPhlan2
summary['Run Assembly']     = params.runAssembly
if (params.runAssembly && params.assembler)   { summary['Assembly tool'] = params.assembler }
summary['Run MetaBAT']      = params.runMetaBAT
// nf-core: Report common parameters here
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch') {
  summary['AWS Region']     = params.awsregion
  summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
  summary['MultiQC maxsize']   = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-hpcmeta-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/hpcmeta Workflow Summary'
    section_href: 'https://github.com/nf-core/hpcmeta'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
// process get_software_versions {
//     publishDir "${params.outdir}/pipeline_info", mode: 'copy',
//         saveAs: { filename ->
//             if (filename.indexOf(".csv") > 0) filename
//             else null
//         }
// 
//     output:
//     file 'software_versions_mqc.yaml' into software_versions_yaml
//     file "software_versions.csv"
// 
//     script:
//     // TODO nf-core: Get all tools to print their version number here
//     """
//     echo $workflow.manifest.version > v_pipeline.txt
//     echo $workflow.nextflow.version > v_nextflow.txt
//     fastqc --version > v_fastqc.txt
//     multiqc --version > v_multiqc.txt
//     scrape_software_versions.py &> software_versions_mqc.yaml
//     """
// }


/*
* Step 1. Trim Reads
*/

process fastqc_raw {
      tag "fastqc_raw_${id}"
      publishDir "${params.outdir}/11-multiqc", mode: 'copy'
  
      input:
      set val(id), file(reads) from read_files_qc
  
      output:
      set val(id), file("*.zip") optional true into RawReadsToMultiQC
    
      """
      fastqc -t ${task.cpus} --noextract ${reads}
      """      
}


 if (!params.skipTrim && params.inputType == 'onp') {
 
     process porechop_ONP {
         tag "porechop_ONP_${id}"
         publishDir "${params.outdir}/1-TrimReadsPorechop", mode: 'copy'
     
         input:
         set val(id), file(reads) from read_files_trimming
     
         output:
         set val(id), file("*.qualtrim.fastq.gz") into trimmedToFilter, trimmedToFakeSanitize, trimmedToQC
         file("*")
         """
         porechop --discard_middle -i ${reads[0]} --threads ${task.cpus} 2> ${id}.porechop.log | \
            seqtk seq -L 1000 | gzip - > ${id}.qualtrim.fastq.gz
         """       
     }
 }
 
if (params.skipTrim && params.inputType == 'onp') {

// need a fake process to connect the channels

    process fakeTrim {
         tag "fakeTrim_${id}"
         publishDir "${params.outdir}/1-FakeTrim", mode: 'copy'
     
         input:
         set val(id), file(reads) from read_files_trimming
     
         output:
         set val(id), file("*.qualtrim.fastq.gz") into trimmedToFilter, trimmedToFakeSanitize, trimmedToQC
         file("*")
         """
         cp ${reads[0]}  ${id}.qualtrim.fastq.gz
         """       
    
    }
}

process fastqc_post {
      tag "fastqc_trimmed_${id}"
      publishDir "${params.outdir}/11-multiqc", mode: 'copy'
  
      input:
      set val(id), file(reads) from trimmedToQC
  
      output:
      set val(id), file("*.zip") optional true into TrimmedReadsToMultiQC
    
      """
      fastqc -t ${task.cpus} --noextract ${reads}
      """      
}

/*
 * Step 2: Remove host DNA
 */
 
// KneadData is not working properly. It produces paired reads that are malformed and cannot be used by bwa
// See post: https://forum.biobakery.org/t/paired-end-data-results-in-unpaired-output/928

if (params.minimapRemoveHost && params.host) {

     process minimap_sanitize_host {
         tag "minimap_sanitize_host_${id}"
         publishDir "${params.outdir}/2-MinimapSanitizeHost", mode: "copy"
 
         input:
         set id, file(reads)  from trimmedToFilter
         file gf from genomeFile
 
         output:
         set id, file("*_sanitized.fastq") into filteredToAssembly,filteredToKraken,filteredToMinimap,filteredToQC,filteredToMetaPhlan2
         file("*sanitized*")

         script: 
         miniparams = params.minimapOpts ? params.minimapOpts : ''
         """
         minimap2 ${miniparams} ${gf} ${reads[0]} | samtools view -b - > ${id}_aln.bam         
         samtools sort -@ ${task.cpus} -o ${id}_aln_sorted.bam  ${id}_aln.bam 
         samtools view -b -f 4 ${id}_aln_sorted.bam > ${id}_sanitized.bam
         samtools fastq -@ ${task.cpus} ${id}_sanitized.bam -n > ${id}_sanitized.fastq

         """
    }
}


if (!params.host && params.inputType == 'onp') {

// need a fake process to connect the channels

    process fakeSanitize {
         tag "fakeSanitize-${id}"
         publishDir "${params.outdir}/2-fakeSanitize", mode: "copy"
 
         input:
         set id, file(reads)  from trimmedToFakeSanitize
 
         output:
         set id, file("*_trimmed.fastq.gz") into filteredToAssembly,filteredToKraken,filteredToMinimap,filteredToQC,filteredToMetaPhlan2
         
         """           
         cp ${reads[0]} ${id}_fsanitized.fastq.gz 
         """   
    }
}

process fastqc_Sanitized {
      tag "fastqc_Sanitized_${id}"
      publishDir "${params.outdir}/11-multiqc", mode: 'copy'
  
      input:
      set val(id), file(reads) from filteredToQC
  
      output:
      set val(id), file("*.zip") optional true into SanitizedReadsToMultiQC
    
      """
      fastqc -t ${task.cpus} --noextract ${reads}
      """      
}


/*
 * Step 4: Phylogenetic profiling
 */


if (params.runMetaPhlan2 && params.inputType == 'onp') {

    process metaphlan2 {
        tag "Metaphlan2-${id}"
        label 'process_high'
        publishDir "${params.outdir}/3-MetaPhlan2", mode: "copy"

        input:
        set val(id), file(freads) from filteredToMetaPhlan2

        output:
        file("${id}/*")

        script:
        runMetaPhlan2params = params.runMetaPhlan2_opts ? params.runMetaPhlan2_opts : ''
        """
        metaphlan2.py ${runMetaPhlan2params} \\
            --bowtie2out ${id}.bowtie2.bz2 \\
            --nproc ${task.cpus} \\
            --input_type fastq \\
            --biom {id}.biom \\
            --tmp_dir /scratch \\
            -o ${id}.profile.txt \\
             ${freads}
        """
    }
}

/*
 * Step 5: taxonomic classification
 */


 if (params.runKraken2 && params.Kraken2DB && params.inputType == 'onp') {

     kraken2DB = file(params.Kraken2DB)

     process kraken2_reads {
         tag "kraken2-${id}"
         publishDir "${params.outdir}/4-Kraken2-Reads", mode: "copy"
 
         input:
         set id, file(freads) from filteredToKraken
 
         output:
         file("${id}*kraken2*")

 
         """
         kraken2 \\
	       --db ${kraken2DB} \\
	       --threads ${task.cpus} \\
	       --output ${id}_reads.kraken2 \\
         --report ${id}_reads.kraken2.report \\
	       ${freads[0]}
         """
     }
 
 }

/*
 * Step 6: Assembly
 */


if (params.runAssembly && params.inputType == 'onp') {

    if (params.assembler == 'megahit') {
    
        process megahit {
            tag "MEGAHIT-${id}"
            label 'process_high'
            publishDir "${params.outdir}/5-MEGAHIT", mode: "copy"
            scratch     '/scratch'
            stageOutMode  'copy'
         
            input:
            set val(id), file(freads) from filteredToAssembly

            output:
            file("${id}/*")
            set val(id), file("*assembly.fasta") into assembly2diamond,assembly2minimap,assembly2MetaBAT2,assembly2Blobtools,assemblyToQC,assemblyToKraken
            set val(id), file("*csv") into assemblyToMultiQC
    
            script:
            megahitparams = params.megahit_opts ? params.megahit_opts : ''
            """
            megahit \\
                -r ${freads[0]} \\
                -t ${task.cpus} \\
                -o ${id} ${megahitparams}
                
            cp ${id}/final.contigs.fa ${id}_assembly.fasta
            
            /home/groups/hpcbio/apps/FAlite/assemblathon_stats.pl -csv ${id}_assembly.fasta           
            """
        }
    }  

    if (params.assembler == 'metaspades') {

        process metaspades {
            tag "metaSPAdes-${id}"
            label 'process_high'
            publishDir "${params.outdir}/5-metaSPADES", mode: "copy"
            scratch     '/scratch'
            stageOutMode  'copy'
         
            input:
            set val(id), file(freads) from filteredToAssembly

            output:
            file("${id}/*")
            set val(id), file("*assembly.fasta") into assembly2diamond,assembly2minimap,assembly2MetaBAT2,assembly2Blobtools,assemblyToQC,assemblyToKraken
            set val(id), file("*csv") into assemblyToMultiQC
                
            script:
            metaspadesparams = params.metaspades_opts ? params.metaspades_opts : ''
            """
            metaspades.py -t ${task.cpus} \\
                --nanopore ${freads[0]}  \\
                -o ${id} ${metaspadesparams}
                
            cp ${id}/scaffolds.fasta ${id}_assembly.fasta
            
            /home/groups/hpcbio/apps/FAlite/assemblathon_stats.pl -csv ${id}_assembly.fasta 
            """
        }
        
    } 
    
    
   if (params.assembler == 'metaflye') {

        process metaflye_ONP {
            tag "metaFlye-${id}"
            publishDir "${params.outdir}/5-metaFlye", mode: "copy"
            scratch     '/scratch'
            stageOutMode  'copy'
    
            input:
            set val(id), file(freads) from filteredToAssembly

            output:
            set val(id), file("*assembly.fasta") into assembly2diamond,assembly2minimap,assembly2MetaBAT2,assembly2Blobtools,assemblyToQC,assemblyToKraken
            set val(id), file("*csv") into assemblyToMultiQC
            file("*")    

            script:
            """
            flye --nano-hq  ${freads[0]} --out-dir ./output --meta --threads ${task.cpus}
            
            cp output/assembly.fasta ${id}_assembly.fasta

            /home/groups/hpcbio/apps/FAlite/assemblathon_stats.pl -csv ${id}_assembly.fasta
            """
        }    
    }
}

/*
 * Step 7: Diamond and or Kraken2 on assembly
 */
 
 if (params.runKraken2 && params.Kraken2DB && params.runAssembly) {

     kraken2DB = file(params.Kraken2DB)

     process kraken2_asm {
         tag "kraken2_asm_${id}"
         publishDir "${params.outdir}/6-Kraken2-Assembly", mode: "copy"
         scratch     '/scratch'
         stageOutMode  'copy'
             
         input:
         set id, file(contigs) from assemblyToKraken
 
         output:
         file("${id}*kraken2*")

 
         """
         kraken2 \\
	       --db ${kraken2DB} \\
	       --threads ${task.cpus} \\
	       --output ${id}_assembly.kraken2 \\
         --report ${id}_assembly.kraken2.report \\
	       $contigs
         """
     }
 }

 
 if (params.runDiamond && params.diamondDB && params.inputType == 'onp') {
    
        diamondDB = file(params.diamondDB)
        // note memory usage; DIAMOND typically requires considerable memory esp. 
        // for larger assemblies and databases (like nr)
        
        process diamond {
            tag "DIAMOND-${id}"
            publishDir "${params.outdir}/7-DIAMOND", mode: "copy"
            scratch     '/scratch'
            stageOutMode  'copy'
                 
            input:
            set val(id), file(contigs) from assembly2diamond

            output:
            set val(id), file("*.daa") into diamondToBlobtools,diamondToView
            file("*.log")
    
            script:
            """
            diamond blastx -p ${task.cpus} \\
                 -d ${diamondDB} \\
                 -q $contigs \\
                 -o ${id}.daa \\
                 --outfmt 100 \\
                 --tmpdir /dev/shm \\
                 --long-reads \\
                 --top 5 \\
                 --evalue 1e-5 \\
                 -v 2> "${id}.log"
            """
        }

        process Diamond_view {
            tag "diamond-view-${id}"
            publishDir "${params.outdir}/7-DIAMOND", mode: "copy"
        
            input:
            set val(id), file(DAA) from diamondToView

            output:
            file("*")
    
            script:
            """
            diamond view -a ${DAA} -f 6 -o ${id}_matches.m8
            """
        }
        
}

/*
 * Step 8: Align reads to asssembly
 */
    
if (params.inputType == 'onp' && params.runMetaBAT || params.runBlobtools) {

        process minimap_aln_asm {
            tag "minimap_aln_asm-${id}"
            publishDir "${params.outdir}/8-minimap_aln_asm", mode: "copy"
        
            input:
            set val(id), file(contigs) from  assembly2minimap         
            set val(id), file(freads) from filteredToMinimap

            output:            
            set val(id), file("${id}*_sorted.bam") into minimap2MetaBAT2,minimap2Blobtools
            set val(id), file("${id}*.bai")  into minimapIdx2MetaBAT2,minimapIdx2Blobtools
            file("${id}*.stats") into alnReads2AsmToQC 
            file("*")      
      
            """
            minimap2 -ax map-ont $contigs ${freads[0]} | samtools view -b - > ${id}_aln_asm.bam       
            samtools sort -@ ${task.cpus} -o ${id}_aln_asm_sorted.bam  ${id}_aln_asm.bam 
            samtools index ${id}_aln_asm_sorted.bam 
            samtools stats ${id}_aln_asm_sorted.bam  > ${id}_aln_asm_sorted.samtools.stats  
            samtools coverage ${id}_aln_asm_sorted.bam > ${id}_aln_asm_sorted.samtools.coverage
            """
        }
}


/*
 * Step 9: Taxonomic binning
 */
 
            
if (params.inputType == 'onp' && params.runMetaBAT) {

        process metabat2_checkM {
            tag "MetaBAT_checkM-${id}"
            publishDir "${params.outdir}/9-MetaBAT-checkM", mode: "copy"
        
            input:
            set val(id), file(aln) from minimap2MetaBAT2
            set val(id2), file(contigs) from assembly2MetaBAT2

            output:
            file("*")
                
            script:
            metabat2params = params.metabat2_opts ? params.metabat2_opts : ''
            """
            echo step one calculate depth
            
            jgi_summarize_bam_contig_depths \\
                --outputDepth ${id}.jgi_depth.txt $aln

            echo step two calculate bins with metabat2

            mkdir  ./${id}_Metabat2_bins/
                     
            metabat2 -i $contigs \\
                -t ${task.cpus} \\
                --unbinned ${metabat2params}  \\
                -a ${id}.jgi_depth.txt \\
                -o ./${id}_Metabat2_bins/

            echo step three run checkm
            
            checkm lineage_wf \\
                -t ${task.cpus} \\
                -x fa \\
                ./${id}_Metabat2_bins/ ./${id}_CheckM > ${id}_CheckM.log

            echo step four generate checkm plots
 
            checkm qa -t 4 -o 2 --tab_table -f ./${id}_CheckM/CheckM_qa_report.tsv ./${id}_CheckM/lineage.ms ./${id}_CheckM           
 
            checkm nx_plot  -x fa ./${id}_Metabat2_bins/ ./${id}_CheckM_nx_plots/
            
            checkm len_hist -x fa ./${id}_Metabat2_bins/ ./${id}_CheckM_lenHist_plots/
                   
            #checkm marker_plot -x fa ./${id}_CheckM ./${id}_Metabat2_bins/ ./${id}_CheckM_marker_plots/

            """
        }
                
}

/*
 * Step 10: Blobtools
 */
 
if (params.runBlobtools) {
    
    process blobtools {
            tag "blobtools-${id}"
            publishDir "${params.outdir}/10-Blobtools", mode: "copy"
            
            input:
            set val(id), file(GENOME) from assembly2Blobtools 
            set val(id), file(DAA) from diamondToBlobtools                    
            set val(id), file(bam) from minimap2Blobtools
            set val(id), file(bai) from minimapIdx2Blobtools
            file PHYLUM  from blobPhylumFile 
            file KINGDOM from blobKingdomFile
            file TAXLIST from blobTaxFile  
                   
            output:
            file("*")
            
            script:
            """ 
            echo step 1 generate TAGC file
            
            perl ${params.blobscript} ${TAXLIST} ${DAA} ${DAA}.tagc
          
                        
            echo step 2 blobtools create
            
            blobtools create -i ${GENOME} -b ${bam} -t ${DAA}.tagc -o ${id}.blobplot
            
            echo step 3 blobtools plot
            
            blobtools plot -i ${id}.blobplot.blobDB.json --colours ${PHYLUM} -p 20

            blobtools plot -i ${id}.blobplot.blobDB.json -r superkingdom --colours ${KINGDOM} -p 20
            
            echo step 4 blobtools view            

            blobtools view -i ${id}.blobplot.blobDB.json -o ${id}.blobplot.blobDB.json_phylum 

            blobtools view -i ${id}.blobplot.blobDB.json -o ${id}.blobplot.blobDB.json_superkingdom -r superkingdom            
                    
            """
    }

}

/*
 * Step x: MultiQC
 */

process multiqc_summary {
    publishDir "${params.outdir}/11-multiqc", mode: 'copy'

    input:
    file ('Raw/*')            from RawReadsToMultiQC.collect()
    file ('Trimmed/*')        from TrimmedReadsToMultiQC.collect()
    file('SanitizedMerged/*') from SanitizedReadsToMultiQC.collect()
    file('assembly/*')        from assemblyToMultiQC.collect()
    file('alnReads2Asm/*')    from alnReads2AsmToQC.collect()

    output:
    file "*"

    script:    
    """
    multiqc . -f -d
    """
}

/*
 * STEP 3 - Output Description HTML
 */
 
// process output_documentation {
//     publishDir "${params.outdir}/pipeline_info", mode: 'copy'
// 
//     input:
//     file output_docs from ch_output_docs
// 
//     output:
//     file "results_description.html"
// 
//     script:
//     """
//     markdown_to_html.r $output_docs results_description.html
//     """
// }

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/hpcmeta] Successful: $workflow.runName"
    if (!workflow.success) {
      subject = "[nf-core/hpcmeta] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/hpcmeta] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/hpcmeta] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if ( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/hpcmeta] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[nf-core/hpcmeta] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if (!output_d.exists()) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[nf-core/hpcmeta]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/hpcmeta]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/hpcmeta v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
