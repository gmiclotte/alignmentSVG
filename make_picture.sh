#!/bin/bash

# usage:
# make_picture.sh <output_name> <contig_fasta> <contig_id> <kmer_count_table> <start_position> <breakpoint_position>

# name for output SVG file
name=${1}
# contig fasta
contig=${2}
# contig index
contig_idx=${3}
# kmer count table
kmer_count=${4}
# start position of the SVG
start=${5}
# position of the breakpoint
breakpoint=${6}

# k-mer size
k=21

# location of the python script
alignmentSVG_py=
# location of the GNUplot file
plotkmers=
# input directory
input_dir=.
# output directory
output_dir=.
# length of the view range, only 100 is supported out of the box
len=100

# tracks
# add as many as you like here
tracks[0]="Uncorrected"
tracks[1]="${input_dir}/mapped.fastq"
tracks[2]="Corrected"
tracks[3]="${input_dir}/mappedCorrected.fastq"
#tracks[4]="Another Corrected"
#tracks[5]="${input_dir}/mappedCorrected.fastq"

kmer_count="${input_dir}/${kmer_count}"
kmer_temp=${output_dir}/${k}_${start}_${len}_results.txt
temp1="${output_dir}/temp1.svg"
temp2="${output_dir}/temp2.svg"

tail -n +$((${start}+1)) "${kmer_count}" | head -n ${len} > "${kmer_temp}"
cut -f1 -d "" "${kmer_temp}" | sort -nk3 | tail -n1
ymax=$((($(cut -f1 -d "" "${kmer_temp}" | sort -nk3 | tail -n1 | awk '{print $3}')/50+1)*50))
keyxpos=$((${start} + 25))
keyypos=$((${ymax} - ${ymax} / 10))

end=$((${start}+${len}))
gnuplot4 -e "kmer_count='${kmer_count}'" -e "ymax='${ymax}'" -e "start='${start}'" -e "end='${end}'" -e "keyxpos='${keyxpos}'" -e "keyypos='${keyypos}'" -e "breakpoint='${breakpoint}'" -e "temp1='${temp1}'" "${plotkmers}"

python "${alignmentSVG_py}" "SAM" "${contig}" "${contig_idx}" "${start}" "${len}" "${kmer_temp}" "${input_dir}/sorted.sam" "${tracks[@]}" > "${temp2}"
str="<g style=\"fill:none; color:red; stroke:currentColor; stroke-width:1.00; stroke-linecap:butt; stroke-linejoin:miter\"><g transform=\"translate(275.5,21.8)\" style=\"stroke:none; fill:black; font-family:Arial; font-size:11.00pt; text-anchor:end\"><text>${k}-mer coverage</text></g><path stroke='rgb( 27, 158, 119)' d='M283.8,17.7 L326.0,17.7'/></g>"
echo $(head -n10 "${temp1}") $(grep -v svg "${temp2}") $str $(tail -n+10 "${temp1}") > "${output_dir}/${name}.svg"
sed -i 's/></>\n</g' "${output_dir}/${name}.svg"
sed -i 's/ width=".*" height="100" viewBox="0 0 .* 100"//g' "${output_dir}/${name}.svg"
svgo "${output_dir}/${name}.svg" "${output_dir}/${name}.min.svg"
xdg-open "${output_dir}/${name}.min.svg" &
