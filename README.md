# Intro

TEnest (https://pubmed.ncbi.nlm.nih.gov/18032588/) is a tool for finding and
annotating transposable element (TE) insertion sites, including "nested" TE
insertions

![](img/1.png) Example of TEnest output

## Organism Repeat Databases:

- [WHEAT](./WHEAT) (from archive.org
  [WHEAT.tar.gz](https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/TE_nest/WHEAT.tar.gz))
- [RICE](./RICE) (from archive.org
  [RICE.tar.gz](https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/TE_nest/RICE.tar.gz))
- [BARLEY](./BARLEY) (from archive.org
  [BARLEY.tar.gz](https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/TE_nest/BARLEY.tar.gz))
- [MAIZE](./MAIZE) (from archive.org
  [source](https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/TE_nest/MAIZE.tar.gz))

# Overview:

TE nest is a software package for annotation and display of nested transposable
element (TE) insertions. TE nest uses a pre-built repeat database for
identification of TEs. TE nest displays the nested TE structure based on an
estimated age since insertion calculated from the divergence of the paired long
terminal repeats (LTRs) of LTR retrotransposons.

- TE nest: Transposable Element Annotation:
  - TE nest (TE nest) searches your input sequence against the chosen repeat
    database to identify nested TE insertions. TE nest will output coordinates
    of the TE annotations, and input these annotations to svg_ltr to provide a
    display showing the chronological order of the nested insertions.
- svg_ltr: Transposable Element Insertion Display:
  - svg_ltr (svg_ltr) gives the ability to separately run the display on TE nest
    annotation files.
- LTR Retrotransposon Insertion Age:
  - TE nest calculates the age since insertion of LTR retrotransposons. This
    calculation is based on the sequence difference between the left and right
    LTRs, when a LTR retrotransposon inserts into the genome the LTRs are
    identical, over evolution the LTRs accumulate mutations. The rate of these
    mutations (rate)(ref) is used to calculate the LTR retrotransposon time
    since insertion (insert calculation).
- TE Database Submission:
  - TE database submission is not yet online. Contact me for more information.

# TE_nest

## Sequence input

Input a DNA sequence into the submission box, or upload a sequence. Your input
sequence can be in fasta format (with a '>header') or simply a plain DNA
sequence. This sequence should be in a flatfile format, with no trailing tabs or
spaces after each line.

## Repeat database

Pick a repeat database from the pull down menu. Alternatively you can choose to
build a custom repeat database.

To use a custom repeat database:

- Select 'custom' from the repeat database pull down menu.
- Upload the three repeat database files. In each case, a multi-sequence fasta
  file is required. This database needs to be 'unique' in that a consensus
  sequence for each repeat family type, or a representitive of each repeat type
  is to be entered. In rare cases using 2-3 entries for a very divergent TE
  family is acceptable. The header of each sequence will become the name of the
  repeat. This repeat name must not contain spaces and case must be the same for
  each time it is entered. For example, an LTR retrotransposon will have three
  sequence entries and one TE list entry; the LTR file, the retrotransposon
  file, the TE file, and the TE list. In each entry the name (header) must be
  the same.
  - LTR retrotransposon sequence file: One sequence entry for each type of LTR
    retrotrasnposon.
  - LTR sequence file: A single LTR of each LTR retrotransposon entered into the
    LTR retrotransposon file.
  - TE sequence file: All TE sequences to be identified. This includes LTR
    retrotransposons, non-LTR retrotransposons, DNA transposons, etc.
- Enter or upload a te_list file

## TE list:

TE nest uses the TE list to find TE annotations in the input sequence.
Customization of the TE list can be used to optimize TE annotations by selecting
or de-selecting certain elements from a provided repeat database. To make a
custom TE list, paste in one of the following organism TE lists, and change to
fit your TE nest run. First column is name of TE, second column is TE type; 0
for LTR retrotransposon, 1 for any other type. LTR retrotransposons (0) must be
at the top of the list.

- Maize TE list
- Rice TE list
- Barley TE list
- Wheat TE list

## TE nest Program Parameters:

TE nest allows customization with the following options:

### Alignment Options:

Alignment options are found three times, once for solo LTRs (SOLOs), once for
LTR retrotransposon middle regions (MIDs), and once for fragmented or whole
non-LTR transposons (FRAG/NLTR). Each of these sections has the same alignment
options, although with different default values. Gap open penalty and gap
extension penalty allow for opening and extending gaps in the alignment.
Alignments reported is how many alignment pieces TE nest will evalulate, raising
this value may find small missing segments, but will greatly increase the
runtime of TE nest. Alignment minimum score is the cutoff value for alignments
to use.

## Output Options and Formats:

Select from the three TE nest output formats:

### Annotation table (\*.LTR file).

The annotation table is the main output of TE nest, giving coordinates for each
identified nested TE. The four annotation types; SOLO (solo LTRs), PAIR (full
LTR retrotransposons), FRAG (fragmented TEs of all types), NLTR (Full length
non-LTR containing TEs). Each the SOLO, FRAG, and NLTR annotation types follow
the same format, one header line and one coordinate line (shown below) per TE
annotation. The PAIR annotation type contains four lines per TE; one header
line, one left LTR coordinate line, one right LTR coordinate line, one middle
region LTR line (shown below).

|                                                             |                                                                                   |
| ----------------------------------------------------------- | --------------------------------------------------------------------------------- |
| Example TE nest \*.LTR annotation file:                     | Explanation:                                                                      |
| 33610                                                       | Sequence length                                                                   |
| input.fasta                                                 | input sequence name                                                               |
| SOLO s0 milt 1 0 20 0                                       | Solo entry, solo number, type, direction, nest group, nest order, nest level      |
| s0 100 276 1 176 662 1140 177 656                           | Solo number, Sequence start, seq end, TE start, TE end, -gap- etc.                |
| SOLO s1 gyma 0 0 32 0                                       |                                                                                   |
| s1 1620 5918 1 4198                                         |                                                                                   |
| PAIR p0 danelle 0 0.030 41 41 0                             | Pair entry, pair number, type, direction, BSR, nest group, nest order, nest level |
| p0 L 8710 13312 1 4602                                      | Pair number, Left LTR, seq start, seq end, TE start, TE end                       |
| p0 R 19507 24778 1 4601                                     | Pair number, Right LTR, seq start, seq end, TE start, TE end                      |
| p0 M 13313 19506 4603 10796                                 | Pair number, Middle region, seq start, seq end, TE start, TE end                  |
| PAIR p1 opie 0 0.023 39 39 0                                |                                                                                   |
| p1 L 6152 6280 1 128 8391 8709 260 577 24107 24778 580 1251 |                                                                                   |
| p1 R 31184 32427 1 1243                                     |                                                                                   |
| p1 M 24779 31183 1266 7670                                  |                                                                                   |
| FRAG f0 dagaf 0 1 1 0                                       | Frag entry, frag number, type, direction, nest group, nest order, nest level      |
| f0 6281 8390 3256 5366                                      | Frag number, Sequence start, seq end, TE start, TE end, -gap- etc.                |
| FRAG f1 grande 0 5 5 0                                      |                                                                                   |
| f1 25203 26653 1 1450 27748 33610 1444 7306                 |                                                                                   |
| NLTR n0 cacta 0 27 21 1                                     | Non-LTR entry, NLTR number, type, direction, nest group, nest order, nest level   |
| n0 26654 27747 1 1094                                       | NLTR number, Sequence start, seq end, TE start, TE end, -gap- etc.                |
| NLTR n1 isb 0 27 17 1                                       |                                                                                   |
| n1 277 661 1 384                                            |                                                                                   |

### SVG TE insertion display (\*.svg file).

The SVG output gives the graphical nesting display. svg_ltr is a separate
program that uses the annotation coordinates file to produce the picture of
nested TEs. This file is scalable vector graphic (SVG) format, viewable in the
latest firefox version. A vector graphic format is used to give resizing
abilities without losing quality. The SVG display is automatically made with a
TE nest sequence input, and can also be re-run by inputting the \*.LTR
annotation file (svg_ltr Display)

### Repeat masked fasta sequence (\*.mask file).

TE nest can also output the original sequence file with all annotated TEs masked
out with Ns. This gives the user the ability to input the masked sequence into
gene prediction programs.

# svg_ltr Display:

## Annotation Input:

Input or upload the _.LTR annotation file, the output of TE nest. You can use
this separate svg_ltr input to rebuild the insertion nesting display, or alter
coordinates on the annotation output and make a new display. You can also use
this function to display TE insertions not identified with TE nest, arrange your
coordinates in the _.LTR annotation format (TE nest output formats).

## svg_ltr Display Subset of Annotations:

A subset of the input sequence can be displayed with svg_ltr. Input the whole
annotation file and enter the start and stop coordinates of interest. The
coordinate window may intersect TE annotations, you can select to either show
these cut insertions as fragments (check the split coordinates box) or to not
display the cut off TEs (deselect the box).

## svg_ltr Program Parameters:

svg_ltr allows customization of the TE nest display by selecting the following
options:

### Display Options:

You can select to display different types of TE nest annotations. Select or
deselect paired LTR Retrotransposons (PAIRs), solo LTRs (SOLOs), whole non-LTR
TEs (NLTRs), or fragmented TEs (FRAGs). At least one check box must be selected.

### White Out Options:

TE nest displays unidentified 'unique' regions found within TE annotations as
white triangles to show the regions does not match the TE sequence. Display of
the white region is usefull for identification of unknown TEs, incomplete TE
annotations in the database, or areas of rearrangements. This option can be
turned on or off.

### Legend Options:

The legend at the bottom of the svg_ltr display can be customized in the
following manner. The default width of TE names in the legend is 5, this can be
shortened or lengthened for different sized TE displays. The scale bar length
can be changed, helpful for longer input sequences. The number of tic marks
within the scale bar can be altered.

### Coordinate Options:

Coordinates of TE annotations can be displayed on the insertion display. Two
types of coordinates are available TE based coordinates show coordinates of the
TE insertions relative to the consensus sequence database TE. Input sequence
based coordinates show the TE insertions relative to the coordinates of the DNA
sequence. Furthermore, each data type can be turned on or off. Select from
PAIRs, SOLOs, FRAGs, or NLTRs to show any of these TE coordinates.

### Insertion Age Options:

Either age since insertion in million years ago (MYA) or the substitution rate
between the left and right LTR can be shown in the box within paired LTR
retrotransposons.

## Footnote

### Background

I (@cmdcolin) downloaded TEnest from the internet archive and am rehosting it on
github for posterity. The plantGDB cgi-bin instance is no longer alive, and I
thought that this tool makes such unique visualizations, that it deserved to be
kept alive

Software link
https://web.archive.org/web/20100507232312/http://www.public.iastate.edu/~imagefpc/Subpages/software.html

Downloaded from
https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/te_nest.html

Original source code
https://web.archive.org/web/20100508181248/http://www.public.iastate.edu/~imagefpc/Subpages/TE_nest/TEnest_scripts.tar.gz

I also looked at two more archive.org links from two later dates but the source
code was unaltered. I believe this is TEnest v1, downloaded from a 2010
archive.org link.

### License

As far as I (@cmdcolin) knows, there is no license attached to the source code

### Concerns

I am just hosting this for software historical purposes. If you wish for this
repository to be taken down, contact me and I will do so!
