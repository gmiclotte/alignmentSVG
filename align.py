'''
	Copyright (C) 2016 Giles Miclotte (giles.miclotte@intec.ugent.be)
	This file is part of alignmentSVG.

	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the
	Free Software Foundation, Inc.,
	59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
'''
from Bio import pairwise2

match = 1
mismatch = -4
gap_open = -6
gap_extend = -1


def extend_cigar(cigar, new):
    count = cigar[-1][0]
    curr = cigar[-1][1]
    if curr == new:
        cigar[-1][0] += 1
    if curr != new:
        cigar += [[1, new]]
    return cigar


def NW(ref, seq):
    alignments = pairwise2.align.globalms(
        ref, seq, match, mismatch, gap_open, gap_extend)
    refal = alignments[0][0]
    seqal = alignments[0][1]
    cigar = [[0, 'M']]
    for i in range(len(refal)):
        r = refal[i]
        s = seqal[i]
        if r == '-':
            cigar = extend_cigar(cigar, 'I')
        elif s == '-':
            cigar = extend_cigar(cigar, 'D')
        else:
            cigar = extend_cigar(cigar, 'M')
    return cigar
