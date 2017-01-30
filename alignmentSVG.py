'''
        Copyright (C) 2016 Giles Miclotte (giles.miclotte@intec.ugent.be)

        This program is free software; you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation; either version 2 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program; if not, write to the
        Free Software Foundation, Inc.,
        59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'''

from __future__ import print_function
import sys
import re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class SAMEntry:
        def __init__(self, line):
                props = line.strip().split('\t')
                self.qname = props[0]
                self.flag = props[1]
                self.rname = props[2]
                self.pos = int(props[3]) - 1 #change to 0-based
                self.mapq = props[4]
                self.cigar = re.findall(r'(\d+)([A-Z]{1})', props[5])
                self.cigar = [(int(c), s) for c, s in self.cigar]
                self.rnext = props[6]
                self.pnext = props[7]
                self.tlen = props[8]
                self.seq = props[9]
                if self.cigar[0][1] == 'S':
                        self.seq = self.seq[self.cigar[0][0]:]
                if self.cigar[-1][1] == 'S':
        	        self.seq = self.seq[:len(self.seq) - self.cigar[-1][0]]
                self.qual = props[10]
                self.optional = {}
                for i in range(11, len(props)):
                        field = props[i].split(':')
                        if field[1] == 'i':
                                field[2] = int(field[2])
                        elif field[1] == 'f':
                                field[2] = float(field[2])
                        self.optional[field[0]] = [field[1], field[2]]
                        
        def len(self):
                return len(self.seq)

class Sequence:
        def __init__(self, meta, seq):
                self.meta = meta
                self.name = meta.strip().split('\t')[0].split(' ')[0]
                self.seq = seq

class XmapEntry:
        def __init__(self):
                self.Alignment = []
        
        def set_props(self, line):
                props = line.strip().split('\t')
                self.XmapEntryID = int(props[0]) - 1
                self.QryContigID = int(props[1]) - 1
                self.RefContigID = int(props[2]) - 1
                self.QryStartPos = float(props[3])
                self.QryEndPos = float(props[4])
                self.RefStartPos = float(props[5])
                self.RefEndPos = float(props[6])
                self.Orientation = props[7] == '+'
                self.Confidence = float(props[8])
                self.HitEnum = props[9]
                self.QryLen = float(props[10])
                self.RefLen = float(props[11])
                self.LabelChannel = props[12]
                for a in props[13][1:-1].split(')('):
                        self.Alignment += [a.split(',')]

def fasta_parse(ifname):
        infile = open(ifname)
        meta = ""
        seq = ""
        for line in infile:
                if line.startswith('>'): #header
                        if meta != "":
                                yield [meta, seq]
                                meta = ""
                                seq = ""
                        meta = line[1:]
                else:
                        seq += line.strip()
        if meta != "":
                yield [meta, seq]
        infile.close()

def fastq_parse(ifname):
        infile = open(ifname)
        meta = ""
        seq = ""
        count = 0
        for line in infile:
                if count == 0:
                        if meta != "":
                                yield [meta, seq]
                                meta = ""
                                seq = ""
                        meta = line[1:]
                elif count == 1:
                        seq = line.strip()
                count = (count + 1) % 4
        if meta != "":
                yield [meta, seq]
        infile.close()


def sam_parse(ifname):
        infile = open(ifname)
        sam = []
        for line in infile:
                if line.startswith('@'): #header
                        continue
                else: #entry
                        sam += [SAMEntry(line)]
        infile.close()
        return sam

def cmap_parse(ifname):
        infile = open(ifname)
        cmap = []
        for line in infile:
                if line.startswith('#'):
                        continue
                else:
                        cmap_entry = line.strip().split('\t')
                        refcontigID = int(cmap_entry[0])
                        while refcontigID > len(cmap):
                                cmap += [[]]
                        cmap[refcontigID - 1] += [float(cmap_entry[5])]
        infile.close()
        return cmap

def xmap_parse(ifname):
        infile = open(ifname)
        xmap = []
        for line in infile:
                if line.startswith('#'):
                        continue
                else:
                        entry = XmapEntry()
                        entry.set_props(line)
                        while entry.RefContigID + 1 > len(xmap):
                                xmap += [[]]
                        xmap[entry.RefContigID] += [entry]
        infile.close()
        return xmap

def bnx_parse(ifname):
        infile = open(ifname)
        bnx = {}
        ID=-1
        for line in infile:
                if line.startswith('#'):
                        continue
                if line.startswith('0'):
                        ID = int(line.strip().split('\t')[1]) - 1
                        continue
                if line.startswith('1'):
                        bnx[ID] = line.strip().split('\t')[1:]
                if line.startswith('Q'):
                        continue
        infile.close()
        return bnx

def entry_before(SVG, ref, entry):
        entry_end = entry.pos + entry.len()
        ref_start = ref.pos
        return entry_end + SVG.min_separation < ref_start

def entry_after(SVG, ref, entry):
        return entry_before(SVG, entry, ref)

def entry_overlap(ref, entry):
        return (not entry_before(ref, entry)) and (not entry_after(ref, entry))

def pile_end_entries(SVG, piles, entry):
        for i in range(len(piles)):
                if entry_after(SVG, piles[i][-1], entry):
                        #append entry
                        piles[i] += [entry]
                        return
        #make new pile
        piles += [[entry]]
        return

def pile_entries(SVG, sam):
        piles = []
        for entry in sam:
                if entry.pos < SVG.begin + SVG.dist and SVG.begin < entry.pos + entry.len():
                        pile_end_entries(SVG, piles, entry)
        return piles

class SVGProperties:
        def __init__(self, begin, dist):
                self.begin = begin
                self.dist = dist
                self.min_separation = 5
                self.height = 3000
#                self.zoom = 50
                self.block_size = 10
                self.nick_width = 4
                self.view_range = [self.begin, self.begin + self.dist]
#                self.min_overlap = 1.0 / 8.0
                self.border_distance = 25
#                self.reference_height = 20
                self.reference_distance = 10
#                self.line_height = 3
                self.line_distance = 1
                self.background_style = "\"fill:white;fill-opacity:1.0;\""
                self.reference_style = "\"fill:#b2df8a;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""
                self.line_style = []
                self.line_style += ["\"fill:#a6cee3;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""]
                self.line_style += ["\"fill:#b2df8a;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""]
                self.nick_style = "\"fill:#33a02c;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""
                self.label_style = "\"fill:#ff0000;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""
                self.label2_style = "\"fill:#0000ff;stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\""
                self.depth = 0

def add_svg_track(SVG, x, y, width, height, style):
        SVG.depth = y + height if y + height > SVG.depth else SVG.depth
        return "\n" + "<rect x=\"" + str(x) + "\" y=\"" + str(y) + "\" width=\"" + str(width) + "\" height=\"" + str(height) + "\" style=" + str(style) + "/>"

def add_svg_track_space(SVG, distance):
        SVG.depth += distance

def draw_seq(SVG, x, y, ref, seq, start, cigar, diff_only = True):
        svg = ""
        ref_idx = 0
        seq_idx = start
        cig_idx = 0
        cig_count = start
        if diff_only:
        	if cigar[cig_idx][1] == 'S' or cigar[cig_idx][1] == 'H':
	        	if seq_idx < cigar[cig_idx][0]:
        			seq_idx = cigar[cig_idx][0]
        			cig_count = cigar[cig_idx][0]
        while seq_idx < len(seq) and ref_idx < len(ref):
                if not diff_only: #draw all the things
                        svg += draw_acgt(SVG, x + ref_idx * SVG.block_size, y, seq[seq_idx])
                        ref_idx += 1
                        seq_idx += 1
                else: #only draw differences
                	if cigar[cig_idx][1] == 'D':
                                svg += draw_acgt(SVG, x + ref_idx * SVG.block_size, y, '-')
                                ref_idx += 1
                        elif cigar[cig_idx][1] == 'I':
                                svg += "<rect x=\"" + str(x + ref_idx * SVG.block_size - 1) + "\" y=\"" + str(y) + "\" width=\"" + str(2) + "\" height=\"" + str(SVG.block_size) + "\" style=" + SVG.label_style + "/>"
                                seq_idx += 1
                        else:
                                if seq[seq_idx] != ref[ref_idx]:
                                        svg += draw_acgt(SVG, x + ref_idx * SVG.block_size, y, seq[seq_idx])
                                ref_idx += 1
                                seq_idx += 1
                        cig_count += 1
                        while cig_count >= cigar[cig_idx][0]:
                                cig_count -= cigar[cig_idx][0]
                                cig_idx += 1
                                if cig_idx == len(cigar):
                                        return svg
                                elif cig_idx == len(cigar) - 1:
                                	if cigar[cig_idx][1] == 'S' or cigar[cig_idx][1] == 'H':
                                		return svg
        return svg

def draw_acgt(SVG, x, y, c):
        return "<text x=\"" + str(x + SVG.block_size / 2) + "\" y=\"" + str(y + SVG.block_size - 1) + "\" text-anchor=\"middle\" alignment-baseline=\"middle\">" + str(c) + "</text>"

def reference_partial_svg(SVG, ref):
        if len(ref.seq) < SVG.view_range[1]:
                SVG.dist = len(ref.seq) - SVG.begin
                SVG.view_range = [SVG.begin, len(ref.seq)]
        x = SVG.border_distance
        y = SVG.border_distance
        w = SVG.dist * SVG.block_size
        h = SVG.block_size
        svg = add_svg_track(SVG, x, y, w, h, SVG.reference_style)
        ref_substr = ref.seq[SVG.view_range[0]:SVG.view_range[1]]
        svg += draw_seq(SVG, x, y, ref_substr, ref_substr, 0, SVG.reference_style, False)
        add_svg_track_space(SVG, SVG.reference_distance)
        return svg

def fasta_partial_svg(SVG, ref, fasta, piles, max_depth, track):
        svg = ""
        for pile in piles:
                y = SVG.depth
                h = SVG.block_size
                for e in pile:
                        start = e.pos - SVG.begin
                        trimmed_start = SVG.border_distance + max(start * SVG.block_size, 0)
                        end = start + e.len()
                        trimmed_end = SVG.border_distance + min(end * SVG.block_size, SVG.dist * SVG.block_size)
                        svg += "\n" #separate entries in SVG by empty line to improve readability
                        x = trimmed_start
                        w = trimmed_end - trimmed_start
                        svg += add_svg_track(SVG, x, y, w , h, SVG.line_style[track%2])
                        ref_substr = ref.seq[SVG.view_range[0] + (start if start > 0 else 0):SVG.view_range[0] + (start if start > 0 else 0) + e.len()]
                        svg += draw_seq(SVG, x, y, ref_substr, fasta[e.qname].seq, -start if start < 0 else 0, e.cigar)
                add_svg_track_space(SVG, SVG.line_distance)
                if SVG.depth + SVG.block_size > max_depth:
                        break
        add_svg_track_space(SVG, SVG.reference_distance)
        return svg

def start_partial_svg(SVG):
        svg = svg = "<svg width=\"" + str(SVG.dist * SVG.block_size + 2 * SVG.border_distance) + "\" height=\"" + str(SVG.height) + "\">"
        svg += "\n" + "<rect x=\"0\" y=\"0\" width=\"" + str(SVG.dist * SVG.border_distance + 2 * SVG.border_distance) + "\" height=\"" + str(SVG.height) + "\" style=" + SVG.background_style + "/>"
        return svg

def end_partial_svg():
        svg = "\n" + "</svg>"
        return svg

def make_svg(SVG, ref, sam, reads):
        svg = reference_partial_svg(SVG, ref)
        ref_depth = SVG.depth
        #total - used - bottom border - remaining borders
        remaining_depth = SVG.height - SVG.depth - SVG.border_distance - (len(reads) - 1) * SVG.reference_distance
        piles = pile_entries(SVG, sam)
        for i in range(len(reads)):
                fasta = reads[i]
                #used + fair portion of the remaining
                max_depth = ref_depth + remaining_depth * (i + 1) / len(reads)
                svg += fasta_partial_svg(SVG, ref, fasta, piles, max_depth, i)
        SVG.height = SVG.depth + SVG.border_distance - SVG.reference_distance
        return start_partial_svg(SVG) + svg + end_partial_svg()

def main(argv=None):
        if argv == None:
                argv = sys.argv
        if len(argv) < 6:
                print('Usage: correctionSVG.py reference.fasta ref_begin_idx ref_distance sorted.sam reads1.fastq <reads2.fastq ...>')
                exit()
        #specify input, can be changed to cl options of course
        reference_file = argv[1]
        SVG = SVGProperties(int(argv[2]), int(argv[3]))
        sam_file = argv[4]
        read_files = [argv[i] for i in range(5, len(argv))] # this can be several files, each will create a new alignment track
        #parse input
        references = [Sequence(s[0], s[1]) for s in fasta_parse(reference_file)]
        sam_entries = sam_parse(sam_file)
        # TODO: we don't actually need to parse the fastq files right now, more efficient to save till moment of drawing that track
        reads = []
        for i in range(len(read_files)):
                fasta = {}
                for s in fastq_parse(read_files[i]):
                        sequence = Sequence(s[0], s[1])
                        fasta[sequence.name] = sequence
                reads += [fasta]
        reference = references[0] #TODO we only look at the first reference, this should probably be specifiec in some other way
        print(make_svg(SVG, reference, sam_entries, reads))

main()
