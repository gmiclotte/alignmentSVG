'''
	Copyright (C) 2016 Giles Miclotte (giles.miclotte@intec.ugent.be)
	This file is part of alignmentSVG.

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
	meta = ''
	seq = ''
	for line in infile:
		if line.startswith('>'): #header
			if len(meta) > 0:
				yield [meta, seq]
				meta = ''
				seq = ''
			meta = line[1:]
		else:
			seq += line.strip()
	if len(meta) > 0:
		yield [meta, seq]
	infile.close()

def fastq_parse(ifname):
	infile = open(ifname)
	meta = ''
	seq = ''
	count = 0
	for line in infile:
		if count == 0:
			if len(meta) > 0:
				yield [meta, seq]
				meta = ''
				seq = ''
			meta = line[1:]
		elif count == 1:
			seq = line.strip()
		count = (count + 1) % 4
	if len(meta) > 0:
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
	return bnx

def entry_before(SVG, ref, entry):
	if SVG.type == 'SAM':
		entry_end = entry.pos + entry.len()
		ref_start = ref.pos
		return entry_end + SVG.min_separation < ref_start
	if SVG.type == 'XMAP':
		entry_end = entry.RefEndPos + entry.QryLen - entry.QryEndPos
		ref_start = -ref.QryStartPos + ref.RefStartPos
		return entry_end + SVG.min_separation < ref_start

def entry_after(SVG, ref, entry):
	return entry_before(SVG, entry, ref)

def pile_end_entries(SVG, piles, entry):
	for i in range(len(piles)):
		if entry_after(SVG, piles[i][-1], entry):
			#append entry
			piles[i] += [entry]
			return
	#make new pile
	piles += [[entry]]
	return

def pile_entries(SVG, entries):
	piles = []
	if SVG.type == 'SAM':
		for entry in entries:
			if entry.pos < SVG.begin + SVG.dist and SVG.begin < entry.pos + entry.len():
				pile_end_entries(SVG, piles, entry)
	elif SVG.type == 'XMAP':
		for entry in entries:
			if entry.Orientation and entry.RefStartPos < SVG.begin + SVG.dist * (1.0 - SVG.min_overlap) and SVG.begin + SVG.dist * SVG.min_overlap < entry.RefEndPos:
				pile_end_entries(SVG, piles, entry)
	return piles

class SVG_properties:
	def __init__(self, data_type, begin, dist):
		# zoom location
		self.begin = begin
		self.dist = dist
		# data type specific settings
		self.type = data_type.upper()
		if self.type == 'SAM':
			self.init_sam()
		elif self.type == 'XMAP':
			self.init_xmap()
		else:
			eprint('Unsupported data type:' + self.type)
			exit()
		# colours
		self.reference_colour = '#b2df8a'
		self.line_colour = ['#a6cee3', '#b2df8a']
		self.nick_colour = '#33a02c'
		self.label_colour = '#ff0000'
		self.label2_colour = '#0000ff'
		self.background_style = '\"fill:white;fill-opacity:1.0;\"'
		self.text_style = '\"writing-mode: bt;text-anchor: middle\"'
		# [current|max] height of the picture
		self.depth = 0
		self.height = 3000
		self.max_subtracks = 200
		# draw text or not
		self.text = False
		# image borders
		base_border = 1
		if self.text:
			self.border = {	'left' : base_border + 6 + 2,
				'lefttext' : base_border + 6,
				'right' : base_border,
				'top' : base_border + 2 * 17 + 2,
				'toptext' : base_border + 17,
				'bottom' : base_border
			}
		else:
			self.border = {	'left' : base_border,
					'lefttext' : base_border,
					'right' : base_border,
					'top' : base_border,
					'toptext' : base_border,
					'bottom' : base_border
			}
		self.border = {	'left' : 83,
				'lefttext' : base_border,
				'right' : 60,
				'top' : 454,
				'toptext' : base_border,
				'bottom' : 10
		}
	
	def init_sam(self):
		self.min_separation = 5
		self.block_size = 8
		self.font_size = str(12.0 / 10.0 * self.block_size)
		self.line_height = self.block_size
		self.view_range = [self.begin, self.begin + self.dist]
		self.track_distance = 10
		self.line_distance = 1
		self.type = 'SAM'
	
	def init_xmap(self):
		self.min_separation = 1000
		self.zoom = 100
		self.nick_width = 2
		self.view_range = [self.begin, self.begin + self.dist]
		self.min_overlap = 1
		self.text_shift = 6
		self.reference_height = 20
		self.track_distance = 10
		self.line_height = 6
		self.line_distance = 1
	
	def track_style(self, colour):
		return '\"fill:' + colour + ';stroke:black;stroke-width:0;fill-opacity:1.0;stroke-opacity:1.0\"'

def add_svg_rect(SVG, x, y, width, height, style):
	SVG.depth = y + height if y + height > SVG.depth else SVG.depth
	return '\n' + '<rect x=\"' + str(x) + '\" y=\"' + str(y) + '\" width=\"' + str(width) + '\" height=\"' + str(height) + '\" style=' + str(style) + '/>'

def add_svg_empty_space(SVG, distance):
	SVG.depth += distance


def draw_seq(SVG, x, y, ref, seq, start, cigar, diff_only = True):
	svg = ''
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
				svg += '<rect x=\"' + str(x + ref_idx * SVG.block_size - 1) + '\" y=\"' + str(y) + '\" width=\"' + str(2) + '\" height=\"' + str(SVG.block_size) + '\" style=' + SVG.track_style(SVG.label_colour) + '/>'
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
	return '<text x=\"' + str(x + SVG.block_size / 2) + '\" y=\"' + str(y + SVG.block_size - 1) + '\" font-size=\"' + SVG.font_size + '\" text-anchor=\"middle\" alignment-baseline=\"middle\">' + str(c) + '</text>'

def add_nicks(SVG, cmap, y, h):
	svg = ''
	for nick in cmap:
		if nick < SVG.view_range[0]:
			continue
		if nick > SVG.view_range[1]:
			break
		x = SVG.border['left'] + (nick - SVG.begin - SVG.zoom / 2) / SVG.zoom
		w = SVG.nick_width
		svg += add_svg_rect(SVG, x, y, w, h, SVG.track_style(SVG.nick_colour))
		#svg += '<g writing-mode=\"tb-rl\" fill=\"black\" font-size=\"4\" transform=\"translate(' + str(x) + ')\"> <text y=\"0px\">' + str(nick) + ' ' + str(x) + '</text></g>'
	return svg

def xmap_partial_svg(SVG, cmap, bnx, piles, max_subtracks, track_names, track):
	eprint('Drawing track: ' + track_names[track] + '.')
	svg = ''
	add_svg_empty_space(SVG, SVG.track_distance)
	subtracks = 0
	for pile in piles:
		subtracks += 1
		y = SVG.depth
		h = SVG.line_height
		for e in pile:
			if not e.Orientation:
				continue
			start = -e.QryStartPos + e.RefStartPos - SVG.begin
			trimmed_start = SVG.border['left'] + max(start / SVG.zoom)
			end = e.RefEndPos - SVG.begin + e.QryLen - e.QryEndPos
			trimmed_end = SVG.border['left'] + min(end / SVG.zoom, SVG.dist / SVG.zoom)
			svg += '\n'
			x = trimmed_start
			w = trimmed_end - trimmed_start
			svg += add_svg_rect(SVG, x, y, w , h, SVG.track_style(SVG.line_colour[track]))
			svg += add_nicks(SVG, [cmap[int(nick)] for nick, label in e.Alignment], y, h)
			'''
			for nick, label in e.Alignment:
				location = (cmap[int(nick) - 1] - SVG.begin - SVG.zoom / 2) / SVG.zoom
				if (SVG.border['left'] < location and location < SVG.dist / SVG.zoom + SVG.border['left']):
					x = SVG.border['left'] + location
					w = SVG.nick_width
					svg += add_svg_rect(SVG, x, y, w , h, SVG.track_style(SVG.nick_colour))
			for label in bnx[e.QryContigID]:
				if e.Orientation:
					location = (float(label) + start - SVG.zoom / 2) / SVG.zoom
					if (SVG.border_distance < location and location < SVG.dist / SVG.zoom + 3 * SVG.border_distance):
						x = 2 * SVG.border_distance + location
						w = SVG.nick_width
						svg += add_svg_rect(SVG, x, y + h / 2, w , h / 2, SVG.track_style(SVG.label_colour))
				else:
					label = float(bnx[e.QryContigID][-1]) - float(label)
					location = (end - float(label) - SVG.zoom / 2) / SVG.zoom
					if (SVG.border_distance < location and location < SVG.dist / SVG.zoom + 3 * SVG.border_distance):
						x = 2 * SVG.border_distance + location
						w = SVG.nick_width
						svg += add_svg_rect(SVG, x, y + h / 2, w , h / 2, SVG.track_style(SVG.label2_colour))
			'''
		add_svg_empty_space(SVG, SVG.line_distance)
		if subtracks == max_subtracks:
			break
	#correct for depth after last subtrack
	add_svg_empty_space(SVG, -SVG.line_distance)   
	if SVG.text:
		#track string
		svg += '\n' + '<g writing-mode=\"tb-rl\" fill=\"black\" font-size=\"8\">'
		x = SVG.border['lefttext']
		y = SVG.depth - (subtracks * (SVG.line_distance + SVG.line_height) - SVG.line_distance) / 2
		svg += '\n' + '<text transform=\"translate(' + str(x) + ', ' + str(y) + ')rotate(270)\" style=' + SVG.text_style + '>' + tracks[track] + '</text>'
		svg += '\n' + '</g>'
	return svg

def fasta_partial_svg(SVG, ref, fasta, piles, max_subtracks, track_names, track):
	eprint('Drawing track: ' + track_names[track] + '.')
	svg = ''
	add_svg_empty_space(SVG, SVG.track_distance)
	subtracks = 0
	for pile in piles:
		subtracks += 1
		y = SVG.depth
		h = SVG.block_size
		for e in pile:
			start = e.pos - SVG.begin
			trimmed_start = SVG.border['left'] + max(start * SVG.block_size, 0)
			end = start + e.len()
			trimmed_end = SVG.border['left'] + min(end * SVG.block_size, SVG.dist * SVG.block_size)
			svg += '\n'
			x = trimmed_start
			w = trimmed_end - trimmed_start
			svg += add_svg_rect(SVG, x, y, w , h, SVG.track_style(SVG.line_colour[track%2]))
			ref_start = SVG.view_range[0] + (start if start > 0 else 0)
			ref_end = ref_start + e.len()
			ref_end = ref_end if ref_end < SVG.view_range[1] else SVG.view_range[1]
			ref_substr = ref.seq[ref_start:ref_end]
			eprint(ref_substr)
			svg += draw_seq(SVG, x, y, ref_substr, fasta[e.qname].seq, -start if start < 0 else 0, e.cigar)
		add_svg_empty_space(SVG, SVG.line_distance)
		if subtracks == max_subtracks:
			break
	#correct for depth after last subtrack
	add_svg_empty_space(SVG, -SVG.line_distance)
	return svg

def reference_partial_svg(SVG, ref):
	if SVG.type == 'SAM':
		return reference_partial_svg_sam(SVG, ref)
	elif SVG.type == 'XMAP':
		return reference_partial_svg_xmap(SVG, ref)

def reference_partial_svg_sam(SVG, ref):
	eprint('Drawing reference track.')
	if len(ref.seq) < SVG.view_range[1]:
		SVG.dist = len(ref.seq) - SVG.begin
		SVG.view_range = [SVG.begin, len(ref.seq)]
	x = SVG.border['left']
	y = SVG.border['top']
	w = SVG.dist * SVG.block_size
	h = SVG.block_size
	svg = add_svg_rect(SVG, x, y, w, h, SVG.track_style(SVG.reference_colour))
	ref_substr = ref.seq[SVG.view_range[0]:SVG.view_range[1]]
	svg += draw_seq(SVG, x, y, ref_substr, ref_substr, 0, SVG.track_style(SVG.reference_colour), False)
	return svg

def reference_partial_svg_xmap(SVG, ref):
	eprint('Drawing reference track.')
	if ref[-1] < SVG.view_range[1]:
		SVG.dist = ref[-1] - SVG.begin
		SVG.view_range = [SVG.begin, ref[-1]]
	if SVG.text:
		svg += '\n' + '<g writing-mode=\"tb-rl\" fill=\"black\" font-size=\"8\">'
		nr = 3
		for i in range(nr):
			#positions
			text = str(int((SVG.begin + i * SVG.dist / (nr - 1)) / 1000)) + 'kb'
			x = SVG.text_shift + SVG.border['left'] + i * (SVG.dist / SVG.zoom - SVG.text_shift) / (nr - 1)
			y = SVG.border['toptext']
			svg += '\n' + '<text transform=\"translate(' + str(x) + ', ' + str(y) + ')rotate(270)\" style=' + SVG.text_style + '>' + text + '</text>'
		svg += '\n' + '</g>'
		#reference string
		svg += '\n' + '<g writing-mode=\"tb-rl\" fill=\"black\" font-size=\"8\">'
		x = SVG.border['lefttext']
		y = SVG.border['top'] + SVG.reference_height / 2
		svg += '\n' + '<text transform=\"translate(' + str(x) + ', ' + str(y) + ')rotate(270)\" style=' +  SVG.text_style + '>Ref</text>'
		svg += '\n' + '</g>'
	x = SVG.border['left']
	y = SVG.border['top']
	w = SVG.dist / SVG.zoom
	h = SVG.reference_height
	svg += add_svg_rect(SVG, x, y, w, h, SVG.track_style(SVG.reference_colour))
	svg += add_nicks(SVG, ref, y, SVG.reference_height)
	return svg

def make_svg(SVG, ref, alnms, tracks, track_names):
	svg = reference_partial_svg(SVG, ref)
	ref_depth = SVG.depth
	# max total - used - bottom border - track borders
	remaining_depth = SVG.height - SVG.depth - SVG.border['bottom'] - len(tracks) * SVG.track_distance
	max_subtracks = (remaining_depth + len(tracks) * SVG.line_distance) / (SVG.line_height + SVG.line_distance)
	max_subtracks = max_subtracks if max_subtracks < SVG.max_subtracks else SVG.max_subtracks
	if SVG.type == 'SAM':
		piles = pile_entries(SVG, alnms)
	for i in range(len(tracks)):
		track = tracks[i]
		if SVG.type == 'SAM':
			svg += fasta_partial_svg(SVG, ref, track, piles, max_subtracks, track_names, i)
		elif SVG.type == 'XMAP':
			piles = pile_entries(SVG, alnms[i])
			svg += xmap_partial_svg(SVG, ref, track, piles, max_subtracks, track_names, i)
	SVG.height = SVG.depth + SVG.border['bottom']
	return start_partial_svg(SVG) + svg + end_partial_svg()

def start_partial_svg(SVG):
	if SVG.type == 'SAM':
		w = SVG.dist * SVG.block_size + SVG.border['left'] + SVG.border['right']
	elif SVG.type == 'XMAP':
		w = SVG.dist / SVG.zoom + SVG.border['left'] + SVG.border['right']
	h = SVG.height
	svg = '<svg width=\"' + str(w) + '\" height=\"' + str(h) + '\">'
	svg += '\n' + '<rect x=\"0\" y=\"0\" width=\"' + str(w) + '\" height=\"' + str(h) + '\" style=' + SVG.background_style + '/>'
	return svg
	
def end_partial_svg():
	svg = '\n' + '</svg>'
	return svg
