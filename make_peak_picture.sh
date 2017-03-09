#!/bin/bash

# usage:
# make_picture.sh <output_name> <contig_fasta> <contig_id> <kmer_count_table> <start_position> <breakpoint_position> <ytic1> <ymax1> <ytic2> <ymin2> <ymax2>

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
# ytic 1
ytic1=${7}
# ymax 1
ymax1=${8}
# ytic2
ytic2=${9}
# ymin 2
ymin2=${10}
# ymax 2
ymax2=${11}

# k-mer size
k=21

# location of the python script
alignmentSVG_py=/home/gmiclott/Documents/Tools/alignmentSVG/alignmentSVG.py
# location of the GNUplot file
plotkmers=/home/gmiclott/Documents/Tools/alignmentSVG/plotkmers
plotkmers_peak="${plotkmers}_peak"
#location of GNUplot 4
gnuplot4=gnuplot4
#location of python 2.7
python=python
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
#tracks[5]="${input_dir}/mappedOtherCorrected.fastq"

#compute some paths and variables
kmer_count="${input_dir}/${kmer_count}"
kmer_temp=${output_dir}/${k}_${start}_${len}_results.txt
temp1="${output_dir}/temp1.svg"
temp2="${output_dir}/temp2.svg"
tail -n +$((${start}+1)) "${kmer_count}" | head -n ${len} > "${kmer_temp}"
cut -f1 -d "" "${kmer_temp}" | sort -nk3 | tail -n1
end=$((${start}+${len}))

#run gnuplot4 to make svg
${gnuplot4} -e "kmer_count='${kmer_count}'" -e "ymax='${ymax}'" -e "start='${start}'" -e "end='${end}'" -e "breakpoint='${breakpoint}'" -e "ytic1='${ytic1}'" -e "ytic2='${ytic2}'" -e "ymax1='${ymax1}'" -e "ymin2='${ymin2}'" -e "ymax2='${ymax2}'" -e "temp1='${temp1}'" "${plotkmers_peak}"

#run python script to make svg
${python} "${alignmentSVG_py}" "SAM" "${contig}" "${contig_idx}" "${start}" "${len}" "${kmer_temp}" "${input_dir}/sorted.sam" "${tracks[@]}" > "${temp2}"

#concat both svg files
str="<g transform=\"translate(0,0)\" id=\"gnuplotkey\"><g style=\"fill:none; color:red; stroke:currentColor; stroke-width:1.00; stroke-linecap:butt; stroke-linejoin:miter\"><g transform=\"translate(275.5,21.8)\" style=\"stroke:none; fill:black; font-family:Arial; font-size:11.00pt; text-anchor:end\"><text>${k}-mer coverage</text></g><path stroke='rgb( 27, 158, 119)' d='M283.8,17.7 L326.0,17.7'/></g></g>"
echo $(head -n-2 "${temp2}") $str $(tail -n+10 "${temp1}") > "${output_dir}/${name}.svg"

#sed magic
viewBoxString=$(grep -o 'viewBox="[0-9.]* [0-9.]* [0-9.]* [0-9.]*"' "${output_dir}/${name}.svg" | grep -o "[0-9. ]*")
read -r -a viewBox <<< "${viewBoxString}"
hshift=35
vshift=5
vsubshift=200
newViewBoxString="${viewBox[0]} ${viewBox[1]} $(echo ${viewBox[2]}+${hshift} | bc) $(echo ${viewBox[3]}+${vshift}+${vsubshift} | bc)"
sed -i "s#${viewBoxString}#${newViewBoxString}#g" "${output_dir}/${name}.svg"
sed -i "s#<g transform=\"translate(0,0)\" id=\"fullSVG\">#<g transform=\"translate(${hshift},${vshift})\" id=\"fullSVG\">#g" "${output_dir}/${name}.svg"
sed -i "s#<g transform=\"translate(0,0)\" id=\"alignmentSVG\">#<g transform=\"translate(0,${vsubshift})\" id=\"alignmentSVG\">#g" "${output_dir}/${name}.svg"
sed -i 's#</svg>#</g></svg>#g' "${output_dir}/${name}.svg"
sed -i 's/></>\n</g' "${output_dir}/${name}.svg"

#open the svg
xdg-open "${output_dir}/${name}.svg" &\
