#!/usr/bin/env nextflow

/*
 *  Copyright (c) 2014-2018, Centre for Genomic Regulation (CRG).
 *  Copyright (c) 2014-2018, Jose Espinosa-Carrasco and the respective authors.
 *
 *  This file is part of Pergola.
 *
 *  Pergola is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Pergola is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Pergola.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Jose Espinosa-Carrasco. CB-CRG. March 2017
 *
 * Script to reproduce Pergola paper figures of CB1 mice experiment
 */

params.recordings     = "$baseDir/small_data/mouse_recordings/"
params.mappings       = "$baseDir/small_data/mappings/b2p.txt"
params.mappings_bed   = "$baseDir/small_data/mappings/bed2pergola.txt"
params.phases         = "$baseDir/small_data/phases/exp_phases.csv"
params.phases_long    = "$baseDir/small_data/phases/exp_phases_whole_exp.csv"
params.mappings_phase = "$baseDir/small_data/mappings/f2g.txt"
params.exp_info       = "$baseDir/small_data/mappings/exp_info.txt"
params.tbl_chromHMM   = "$baseDir/small_data/chromHMM/cellmarkfiletable"
params.output         = "files/"
params.image_format   = "tiff"

log.info "Mouse - Pergola - Reproduce  -  version 0.2.0"
log.info "====================================="
log.info "mice recordings          : ${params.recordings}"
log.info "mappings                 : ${params.mappings}"
log.info "mappings bed             : ${params.mappings_bed}"
log.info "experimental phases      : ${params.phases}"
log.info "experimental phases long : ${params.phases_long}"
log.info "mappings phases          : ${params.mappings_phase}"
log.info "experimental info        : ${params.exp_info}"
log.info "output                   : ${params.output}"
log.info "image format             : ${params.image_format}"
log.info "chromHMM config table    : ${params.tbl_chromHMM}"
log.info "\n"

// Example command to run the script
/*
nextflow run mouse-pergola-reproduce.nf \
  --recordings='small_data/mouse_recordings/' \
  --mappings='small_data/mappings/b2p.txt' \
  --mappings_bed='small_data/mappings/bed2pergola.txt' \
  --phases='small_data/phases/exp_phases.csv' \
  --phases_long='small_data/phases/exp_phases_whole_exp.csv' \
  --mappings_phase='small_data/mappings/f2g.txt' \
  --exp_info='small_data/mappings/exp_info.txt' \
  --image_format='tiff' \
  --tbl_chromHMM="small_data/chromHMM_files/cellmarkfiletable" \
  -with-docker
*/

/*
 * Input parameters validation
 */
mapping_file = file(params.mappings)
mapping_bed_file = file(params.mappings_bed)
mapping_file_bG = file(params.mappings)
mapping_file_phase = file(params.mappings_phase)

exp_phases = file(params.phases)
exp_phases_long = file(params.phases_long)
exp_info = file(params.exp_info)

cell_mark_file_tbl = file(params.tbl_chromHMM)

//cell_mark_file_tbl.println()

/*
 * Input files validation
 */
if( !mapping_file.exists() ) exit 1, "Missing mapping file: ${mapping_file}"
if( !mapping_file_phase.exists() ) exit 1, "Missing mapping phases file: ${mapping_file_phase}"
if( !exp_phases.exists() ) exit 1, "Missing phases file: ${exp_phases}"
if( !exp_phases_long.exists() ) exit 1, "Missing long phases file: ${exp_phases_long}"
if( !exp_info.exists() ) exit 1, "Missing experimental info file: ${exp_info}"
if( !cell_mark_file_tbl.exists() ) exit 1, "Missing configuration file for chromHMM: ${cell_mark_file_tbl}"

/*
 * Read image format
 */
image_format = "${params.image_format}"

/*
 * Create a channel for mice recordings
 */
Channel
    .fromPath( "${params.recordings}intake*.csv" )
    .ifEmpty { error "Cannot find any CSV file with mice data" }
    .set { mice_files }

mice_files.into { mice_files_bed; mice_files_bedGraph }

/*
 * Create a channel for mice recordings
 */
Channel
    .fromPath( params.recordings )
    .set { mice_files_preference }

/*
 * Calculates mice preference statistics
 */
process behavior_by_week {
    publishDir "results_for_test/", mode: 'copy', overwrite: 'true'

  	input:
  	file file_preferences from mice_files_preference
  	file mapping_file
    file mapping_file_phase
    file exp_phases
    file exp_phases_long

  	output:
  	file 'behaviors_by_week' into d_behaviors_by_week
  	stdout into max_time
    file 'exp_phases' into exp_phases_bed_to_wr, exp_phases_bed_to_wr2
    file 'exp_phases_sushi' into exp_phases_bed_sushi, exp_phases_bed_gviz

	file 'stats_by_phase/phases_dark.bed' into exp_circadian_phases_sushi, exp_circadian_phases_gviz, days_bed_igv, days_bed_shiny, days_bed_deepTools
    file 'Habituation_dark_long.bed' into bed_dark_habituation
    file 'Development_dark_long.bed' into bed_dark_development

  	"""
  	mice_stats_by_week.py -f "${file_preferences}"/intake*.csv -m ${mapping_file} -s "sum" -b feeding -p ${exp_phases} \
  	                      -pl ${exp_phases_long} -mp ${mapping_file_phase}

  	mkdir stats_by_phase
  	mkdir behaviors_by_week
  	cp exp_phases.bed exp_phases
  	tail +2 exp_phases > exp_phases_sushi
  	mv *.bed  stats_by_phase/
  	mv stats_by_phase/Development_dark_long.bed ./
  	mv stats_by_phase/Habituation_dark_long.bed ./
  	mv *.tbl  behaviors_by_week/
  	"""
}

def igv_files_by_group ( file ) {

    def map_id_group = [ "ctrl" : [1,3,5,7,11,13,15,17],
                         "hf" : [2,4,8,10,12,14,16,18] ]

    def id = file.split("\\_")[1]

    def food = file.split("\\_")[4].split("\\.")[0]

    def ext = file.split("\\.")[1]

    if ( map_id_group.get("ctrl").contains(id.toInteger()) && food == "sc" )
    return "igv/1_ctrl_food_sc/${id}.${ext}"

    if ( map_id_group.get("hf").contains(id.toInteger()) && food == "fat" )
    return  "igv/2_hf_food_fat/${id}.${ext}"

}

/*
 * Converts input behavioral trajectory of mice into BED files (discrete intervals)
 */
process convert_bed {
    publishDir params.output_res, mode: 'copy', pattern: "tr*food*.bed", saveAs: this.&igv_files_by_group
    publishDir "results_bed/", mode: 'copy', pattern: "tr_*.bed", overwrite: 'true'

  	input:
  	file ('batch') from mice_files_bed
  	file mapping_file
  	file mapping_bed_file

  	output:
  	file 'tr*food*.bed' into bed_out, bed_out_shiny_p, bed_out_gviz, bed_out_sushi
  	file 'dir_bed' into dir_bed_to_bin
  	file 'dir_bed_hab' into dir_bed_i_habituation
  	// file 'tr*{water,sac}*.bed' into bed_out_drink
  	file 'phases_light.bed' into phases_night
  	file '*.fa' into out_fasta

  	"""
  	pergola_rules.py -i ${batch} -m ${mapping_file} -f bed -nt -e -bl -dl food_sc food_fat -d all

    shopt -s nullglob

	## ctrl
  	for f in {1,3,,5,7,11,13,15,17}
  	do
  	    mkdir -p work_dir

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
  	        cd work_dir
  	        track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
  	        echo -e "food_sc\tblack" > dict_color
  	        echo -e "food_fat\tblack" >> dict_color
  	        pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
  	        in_f_sc=`ls tr_chr1*food_sc.bed`
            mv "\$in_f_sc" "`echo \$in_f_sc | sed s/chr1/\${f}/`"
  	        cd ..
  	        echo \$PWD
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

    for f in {2,4,8,10,12,14,16,18}
  	do
  	    mkdir -p work_dir

  	    files=( tr_"\$f"_* )

  	    if (( \${#files[@]} )); then
  	        cd work_dir
  	        track_int=`ls ../"tr_"\$f"_"*`
  	        mv \${track_int} track_int
  	        echo -e "food_sc\torange" > dict_color
  	        echo -e "food_fat\torange" >> dict_color
  	        pergola_rules.py -i track_int -m ../${mapping_bed_file} -c dict_color -f bed -nt -e -nh -s 'chrm' 'start' 'end' 'nature' 'value' 'strain' 'color' -dl food_sc food_fat -d all
  	        in_f_sc=`ls tr_chr1*food_sc.bed`
            mv "\$in_f_sc" "`echo \$in_f_sc | sed s/chr1/\${f}/`"
  	        cd ..
  	        mv work_dir/tr*.bed ./
  	        mv work_dir/*.fa ./
        fi
  	done

  	mkdir dir_bed
  	mkdir dir_bed_hab

  	cp tr*.bed ./dir_bed
  	cp tr*.bed ./dir_bed_hab
  	"""
}

hab_bed = Channel.fromPath ("/Users/jespinosa/git/mouse_chrom_hmm/small_data/phases/habituation.bed")

process intersect_bed_habituation {
    publishDir "results_habituation_binning/", mode: 'copy', overwrite: 'true'

    input:
    file dir_bed_data from dir_bed_i_habituation
    file habituation from hab_bed

    output:
    file 'dir_bed_hab/*.bed' into bed_habituation

    """
    mv ${habituation} ${dir_bed_data}/
    cd ${dir_bed_data}

    rm -f habituation.tr*.bed

    for file_bed in tr_*.bed
    do
	    bedtools intersect -a \${file_bed} -b ${habituation} > "habituation."\${file_bed}
    done
    """
}

/*
 * Converts input behavioral trajectory of mice into bedGraph files showing a continuous score along time windows (30 min)
 */
process convert_bedGraph {

    publishDir params.output_res, mode: 'copy', pattern: "tr*food*.bedGraph", saveAs: this.&igv_files_by_group

  	input:
  	file ('batch_bg') from mice_files_bedGraph
  	file mapping_file_bG
  	val max from max_time.first()

  	output:
  	file 'tr*food*.bedGraph' into bedGraph_out, bedGraph_out_shiny_p, bedGraph_out_gviz, bedGraph_out_sushi, bedGraph_out_bigwig
  	file 'chrom.sizes' into chrom_sizes, chrom_sizes_chromHMM_b, chrom_sizes_chromHMM_l
  	//file 'tr*{water,sac}*.bedGraph' into bedGraph_out_drink

  	"""
  	pergola_rules.py -i ${batch_bg} -m ${mapping_file_bG} -max ${max} -f bedGraph -w 1800 -nt -e -dl food_sc food_fat -d all
  	"""
}

/*
 * Binning of meals based on meal duration
 */
n_bins = 5
process bin {

    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file (dir_bed_feeding) from dir_bed_to_bin
    file 'cellmarkfiletable' from cell_mark_file_tbl

    output:
    file 'bed_binned' into dir_bed_binned
    file '*.binned' into cell_mark_file_tbl_binned
    file "meal_length_distro_binned.${image_format}" into plot_distro_binned

    """
    distro_meals_to_bin.R --path_bed_files=${dir_bed_feeding} \
                          --n_bins=${n_bins} \
                          --image_format=${image_format} > bins.txt

    for file_bed in ${dir_bed_feeding}/*.bed
    do
        # bin_length_by_sliding_win.py -b \${file_bed} -ct 1 604800 -bins "\$(< bins.txt)"
        bin_length_by_sliding_win.py -b \${file_bed} -bins "\$(< bins.txt)"
    done

    mkdir bed_binned
    mv *.bed bed_binned

    bins_string="\$(tr -d "\n\r" < bins.txt)"
    IFS=' ' read -r -a bins_ary <<< \$bins_string

    length_ary=\${#bins_ary[@]}
    i_last=\$((length_ary-1))
    i_for=\$((length_ary-2))

    awk -v bin_0=\${bins_ary[0]} \
        '{print \$1"\t0_"bin_0"_"\$2"\t0_"bin_0"_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"

    for index in \$(seq 0 \$i_for); do
        next_i=\$((index+1))
        awk -v bin_1=\${bins_ary[index]} -v bin_2=\${bins_ary[next_i]} \
            '{print \$1"\t"bin_1"_"bin_2"_"\$2"\t"bin_1"_"bin_2"_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    done

    awk -v  bin_l=\${bins_ary[i_last]} \
        '{print \$1"\t"bin_l"_"\$2"\t"bin_l"_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    """
}

/*
 * chromHMM binarizes feeding bed files
 */
process binarize {
    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_b
    file dir_bed_binned from dir_bed_binned
    file 'cellmarkfiletable_binned' from cell_mark_file_tbl_binned

    output:
    file 'output_dir' into output_dir_binarized

    """
    mkdir output_dir

    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_binned} ${cellmarkfiletable_binned} output_dir
    """
}

/*
 * chromHMM learn model
 * In this case we use to learn the model all the data
 */
n_states = 3
process HMM_model_learn {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l
    file 'input_binarized' from output_dir_binarized

    output:
    file 'output_learn/*dense*.bed' into HMM_model_ANNOTATED_STATES
    file 'output_learn/*.*' into HMM_full_results
    file '*.bed' into segmentation_bed

    """
    mkdir output_learn

    # blue 1 active
    echo -e "1\t0,0,255" > colormappingfile
    # red 2 resting
    echo -e "2\t255,0,0" >> colormappingfile
    # yellow 3 snacking
    echo -e "3\t255,255,0" >> colormappingfile
    # green 4
    echo -e "4\t0,255,0" >> colormappingfile

    # java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes}  -printstatebyline test_feeding/output/outputdir test_feeding/output/outputdir_learn ${n_states} test_feeding/input/chrom.sizes
    java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 \
                                                         -l ${chrom_sizes} \
                                                         input_binarized output_learn \
                                                         ${n_states} ${chrom_sizes}

    for dense_file in output_learn/*segments*.bed
    do
        filename=\$(basename -- "\$dense_file")
        filename="\${filename%.*}"
        mice_id=\$(echo \$filename | cut -f2 -d_)

        java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \${dense_file} \${mice_id} \${filename}
    done
    """
}


/*

rm -f [0-9]*.bed

    for file_bed in ${dir_bed_feeding}/*.bed
    do
        # bin_length_by_sliding_win.py -b \${file_bed} -ct 1 604800 -bins 30 120
        bin_length_by_sliding_win.py -b \${file_bed} -ct ${start} ${end} -bins ${bins}
    done
    echo ${index}
    mkdir output_bed_binned
    mv *.bed output_bed_binned
*/

/*
 * chromHMM binarizes feeding bed files
 */
/*
process bin_binarize {
    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_b
    // file dir_bed_feeding from dir_bed_to_chromHMM
    file dir_bed_feeding from dir_bed_to_bin
    file 'cellmarkfiletable' from cell_mark_file_tbl

    output:
    file 'output_dir' into output_dir_binarized

    """
    cd ${dir_bed_feeding}

    rm -f [0-9]*.bed

    for file_bed in *.bed
    do
	    cat \${file_bed} | awk -F'\t' '\$3-\$2<=30 {print \$0}' > "30_\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>30 && \$3-\$2<=120{print \$0}' > "30_120\${file_bed}"
	    cat \${file_bed} | awk -F'\t' '\$3-\$2>120 {print \$0}' > "120_\${file_bed}"
    done

    cd ..

    awk '{print \$1"\t""30_"\$2"\t""30_"\$3}' cellmarkfiletable > "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""30_120"\$2"\t""30_120"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"
    awk '{print \$1"\t""120_"\$2"\t""120_"\$3}' cellmarkfiletable >> "${cellmarkfiletable}.binned"

    mkdir output_dir

    java -mx4000M -jar /ChromHMM/ChromHMM.jar BinarizeBed -b 300 -peaks ${chrom_sizes} ${dir_bed_feeding} "${cellmarkfiletable}.binned" output_dir
    """
}
*/

/*
 * chromHMM learn model
 * In this case we use to learn the model all the data
 */

n_states = 3
process HMM_model_learn {

    publishDir "${params.output_res}/chromHMM", mode: 'copy', overwrite: 'true'

    input:
    file chrom_sizes from chrom_sizes_chromHMM_l
    file 'input_binarized' from output_dir_binarized

    output:
    file 'output_learn/*dense*.bed' into HMM_model_ANNOTATED_STATES
    file 'output_learn/*.*' into HMM_full_results
    file '*.bed' into segmentation_bed

    """
    mkdir output_learn

    # blue 1 active
    echo -e "1\t0,0,255" > colormappingfile
    # red 2 resting
    echo -e "2\t255,0,0" >> colormappingfile
    # yellow 3 snacking
    echo -e "3\t255,255,0" >> colormappingfile

    # java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 -l  ${chrom_sizes}  -printstatebyline test_feeding/output/outputdir test_feeding/output/outputdir_learn ${n_states} test_feeding/input/chrom.sizes
    java -mx4000M -jar /ChromHMM/ChromHMM.jar LearnModel -b 300 \
                                                         -l ${chrom_sizes} \
                                                         input_binarized output_learn \
                                                         ${n_states} ${chrom_sizes}

    for dense_file in output_learn/*segments*.bed
    do
        filename=\$(basename -- "\$dense_file")
        filename="\${filename%.*}"
        mice_id=\$(echo \$filename | cut -f2 -d_)

        java -mx4000M -jar /ChromHMM/ChromHMM.jar MakeBrowserFiles -c colormappingfile \${dense_file} \${mice_id} \${filename}
    done
    """
}

/*
 *
 */
process plot_HMM_states {
    publishDir "results/", mode: 'copy', overwrite: 'true'

    input:
    file 'output_learn/*' from segmentation_bed.collect()

    output:
    file "segmentation_HMM.${image_format}" into plot_HMM_segmentation

    """
    plot_HMM_segmentation.R --path_bed_files=output_learn \
                            --ini_time=0 \
                            # --end_time=1814400 \
                            --image_format=${image_format}
    """
}
