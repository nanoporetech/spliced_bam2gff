package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/feat/gene"
	"github.com/biogo/biogo/feat/genome"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/hts/sam"
)

type Locus struct {
	Chrom string
	Start int
	End   int
	Feats []gff.Feature
	Order int
	Size  int
}

// NewLocus return a pointer to a new Locus structure.
func NewLocus(chrom string, start, end, order int) *Locus {
	return &Locus{chrom, start, end, []gff.Feature{}, order, 0}
}

// String representation of a locus object.
func (l *Locus) String() string {
	return fmt.Sprintf("%09d_%s:%d:%d", l.Order, l.Chrom, l.Start, l.End)
}

// SplicedBam2PartGFF converts spliced BAM alignments to locus-partitioned GFF files.
func SplicedBam2PartGFF(inBam string, outDir string, minBundle int, nrProcBam int, minimapInput bool, strandBehaviour int, maxDel int, keepS bool) {
	err := os.MkdirAll(outDir, 0750)
	if err != nil {
		L.Fatalf("Could not create GFF output directory: %s", err)
	}
	var locus *Locus
	locusCache := make([]*Locus, 0, 1000)
	var count, bundleCount, bc int
	locusChan := SplicedBam2Loci(inBam, nrProcBam, minimapInput, strandBehaviour, maxDel, keepS)
	for locus = range locusChan {
		count += locus.Size
		bc += locus.Size
		if locus.Size > minBundle {
			locusCache = append(locusCache, locus)
		}
		if (bc > minBundle) || (locus.Size > minBundle) {
			bundleName := fmt.Sprintf("%09d_%s:%d:%d_bundle.gff", bundleCount, locus.Chrom, locusCache[0].Start, locus.End)
			outName := path.Join(outDir, bundleName)
			outFh, err := os.Create(outName)
			if err != nil {
				L.Fatalf("Could not create GFF output file: %s", err)
			}
			outBuff := bufio.NewWriter(outFh)
			gffWriter := gff.NewWriter(outBuff, 1000, true)
			for _, l := range locusCache {
				WriteFeatures(l.Feats, gffWriter)
			}

			outBuff.Flush()
			outFh.Close()
			locusCache = locusCache[:0]
			bc = 0
			bundleCount++
		} else {
			locusCache = append(locusCache, locus)
		}
	}
	if len(locusCache) > 0 {
		bundleName := fmt.Sprintf("%09d_%s:%d:%d_bundle.gff", bundleCount, locus.Chrom, locusCache[0].Start, locus.End)
		outName := path.Join(outDir, bundleName)
		outFh, err := os.Create(outName)
		if err != nil {
			L.Fatalf("Could not create GFF output file: %s", err)
		}
		outBuff := bufio.NewWriter(outFh)
		gffWriter := gff.NewWriter(outBuff, 1000, true)
		for _, l := range locusCache {
			WriteFeatures(l.Feats, gffWriter)
		}

		outBuff.Flush()
		outFh.Close()
		L.Println(len(locusCache))
		locusCache = locusCache[:0]
		bc = 0
		bundleCount++

	}
	L.Printf("Written %d transcripts to %d loci and %d bundles.", count, locus.Order+1, bundleCount)
}

// SplicedBam2PartGFF converts spliced BAM alignments to GFF features grouped by loci.
func SplicedBam2Loci(inBam string, nrProcBam int, minimapInput bool, strandBehaviour int, maxDel int, keepS bool) chan *Locus {
	bamReader := NewBamReader(inBam, nrProcBam)
	outChan := make(chan *Locus, 100)
	var cache []*sam.Record
	var loc *Locus
	var chrom string
	locCount := 0

	go func() {
		// Ierate over BAM records:
		for {
			record, err := bamReader.Read()

			if err == io.EOF {
				break
			}

			// Turn mapped SAM records into GFF:
			if record.Flags&sam.Unmapped == 0 {
				if !keepS {
					if (record.Flags&sam.Secondary != 0) || (record.Flags&sam.Supplementary != 0) || hasSupp(record) {
						continue
					}
				}
				if cache == nil {
					cache = make([]*sam.Record, 0, 5000)
					loc = NewLocus(record.Ref.Name(), record.Start(), record.End(), locCount)
					cache = append(cache, record)
					loc.Size++
					locCount++
				} else if (record.Start() > loc.End) || (record.Ref.Name() != loc.Chrom) {
					for _, r := range cache {
						loc.Feats = append(loc.Feats, SplicedSAM2GFF(r, minimapInput, strandBehaviour, maxDel)...)
					}
					outChan <- loc
					loc = NewLocus(record.Ref.Name(), record.Start(), record.End(), locCount)
					if loc.Chrom != chrom {
						chrom = loc.Chrom
						L.Printf("Processing chromosome: %s\n", chrom)
					}
					cache = cache[:0]
					cache = append(cache, record)
					loc.Size++
					locCount++
				} else {
					if record.Start() < cache[len(cache)-1].Start() {
						L.Fatalf("BAM file is not sorted! Offending records: %s %s\n", cache[len(cache)-1].Ref.Name(), record.Ref.Name())
					}
					if record.End() > loc.End {
						loc.End = record.End()
					}
					cache = append(cache, record)
					loc.Size++

				}
			}
		}

		for _, r := range cache {
			loc.Feats = append(loc.Feats, SplicedSAM2GFF(r, minimapInput, strandBehaviour, maxDel)...)
		}
		outChan <- loc
		close(outChan)

	}()
	return outChan
}

// Turn a BAM file containing sliced alignments into GFF2 format annotation.
func SplicedBam2GFF(inBam string, out io.Writer, nrProcBam int, minimapInput bool, strandBehaviour int, maxDel int, keepS bool) {
	bamReader := NewBamReader(inBam, nrProcBam)

	gffWriter := gff.NewWriter(out, 1000, true)
	count := 0

	// Ierate over BAM records:
	for {
		record, err := bamReader.Read()

		if err == io.EOF {
			break
		}

		// Turn mapped SAM records into GFF:
		if record.Flags&sam.Unmapped == 0 {
			if !keepS {
				if (record.Flags&sam.Secondary != 0) || (record.Flags&sam.Supplementary != 0) || hasSupp(record) {
					continue
				}
			}
			SplicedSAM2GFFWrite(record, gffWriter, minimapInput, strandBehaviour, maxDel)
			count++
		}
	}
	L.Printf("Written %d transcripts.", count)
}

// Create a new gene.CodingTranscript object from SAM reference, position and orientation.
func NewCodingTranscript(chrom *sam.Reference, id string, pos int, strand feat.Orientation) *gene.CodingTranscript {

	// This will allocate a new chromosome for each transcript
	// but withing this application that should be OK:
	ch := &genome.Chromosome{
		Chr:      chrom.Name(),
		Desc:     chrom.Name(),
		Length:   chrom.Len(),
		Features: nil,
	}

	tr := &gene.CodingTranscript{
		ID:       id,
		Loc:      ch,
		Offset:   pos,
		Orient:   strand,
		Desc:     id,
		CDSstart: 0,
		CDSend:   0,
	}

	return tr
}

func hasSupp(rec *sam.Record) bool {
	_, ok := rec.Tag([]byte("SA"))
	if ok {
		return true
	}
	return false
}

// Get orientation from transcript strand tag (either XS, or ts for minimap2).
func getTrStrand(rec *sam.Record, minimapInput bool) feat.Orientation {
	var aux sam.Aux
	if minimapInput {
		aux, _ = rec.Tag([]byte("ts"))
	} else {
		aux, _ = rec.Tag([]byte("XS"))
	}

	// We got the tag value:
	if aux != nil {
		// Convert tag value to string:
		strand := string(aux.Value().(uint8))
		// Decide orientation:
		switch strand {
		case "+":
			return feat.Forward
		case "-":
			return feat.Reverse
		case "?":
			return feat.NotOriented
		default:
			L.Fatalf("Unknown orientation string: %s\n", strand)
		}
	} else {
		//L.Printf("Missing strand tag in record: %s\n", rec.Name)
	}

	// Missing tag, feature not oriented:
	return feat.NotOriented
}

// Flip orientation:
func flipOrientation(orient feat.Orientation) feat.Orientation {
	switch orient {
	case feat.Forward:
		return feat.Reverse
	case feat.Reverse:
		return feat.Forward
	case feat.NotOriented:
		return feat.NotOriented
	default:
		L.Fatalf("Unknown orientation: %s", orient)
	}
	return feat.NotOriented
}

// Decide on the feature strand depending on the transcript strand tag and read orientation:
func figureStrand(readStrand, trStrand feat.Orientation, minimapInput bool, strandBehaviour int) feat.Orientation {

	// Use read orientation as feature strand:
	if strandBehaviour == StrandRead {
		return readStrand
	}

	var strand feat.Orientation

	// Strand tag is missing:
	if trStrand == feat.NotOriented {
		switch strandBehaviour {
		case StrandTag:
			strand = feat.NotOriented // Strand tag takes precedence, feature is not oriented.
		case StrandTagRead:
			strand = readStrand // Fallback to read orientation.
		}
		return strand
	}

	// Transript strand tag is present:

	switch minimapInput {
	case true:
		if trStrand == feat.Reverse {
			strand = flipOrientation(readStrand) // Flip orientaton.
		} else {
			strand = readStrand // Use read strand.
		}
	case false:
		// Input is not minimap2, use transcript strand tag as feature orientation.
		strand = trStrand
	}

	return strand
}

// Convert SAM record into GFF2 records. Each read will be represented as a distinct transcript.
func SplicedSAM2GFFWrite(record *sam.Record, gffWriter *gff.Writer, minimapInput bool, strandBehaviour int, maxDel int) {
	WriteFeatures(SplicedSAM2GFF(record, minimapInput, strandBehaviour, maxDel), gffWriter)
}

// Convert SAM record into GFF2 records. Each read will be represented as a distinct transcript.
func SplicedSAM2GFF(record *sam.Record, minimapInput bool, strandBehaviour int, maxDel int) []gff.Feature {

	//Get read strand:
	var readStrand feat.Orientation = feat.Forward
	if record.Flags&sam.Reverse != 0 {
		readStrand = feat.Reverse
	}

	// Get transcript strand:
	trStrand := getTrStrand(record, minimapInput)

	// Decide feature strand:
	strand := figureStrand(readStrand, trStrand, minimapInput, strandBehaviour)

	transcript := NewCodingTranscript(record.Ref, record.Name, record.Pos, strand)

	exons := make(gene.Exons, 0) // To accumulate exons.

	// First exon starts at record position:
	var currBlockStart int = record.Pos
	var currBlockLen int = 0
	var exonNr int = 0

CIGAR_LOOP: // Iterate over CIGAR:
	for _, cigar := range record.Cigar {
		op := cigar.Type()
		length := cigar.Len()

		switch op {

		// Soft clip, hard clip, or insertion - do not consume reference:
		case sam.CigarSoftClipped, sam.CigarHardClipped, sam.CigarInsertion:
			continue CIGAR_LOOP

			// Match, mismatch or deletion - add to current exon length:
		case sam.CigarDeletion:
			if length < maxDel {
				currBlockLen += length
			} else {
				exonStart := currBlockStart              // Previous exon starting here.
				exonEnd := currBlockStart + currBlockLen // Previous exon ends here.

				// Create exon object:
				exonId := fmt.Sprintf("exon_%d", exonNr)
				exon := gene.Exon{transcript, exonStart - record.Pos, exonEnd - exonStart, exonId}

				// Discard zero length exons - FIXME: maybe this should not happen.
				if exon.Len() > 0 {
					exons = append(exons, exon) // Register exon.
				}

				currBlockLen = 0                  // Reset exon length counter.
				currBlockStart = exonEnd + length // Next exon starts after the N operation.
				exonNr++

			}
		case sam.CigarMatch, sam.CigarEqual, sam.CigarMismatch:
			currBlockLen += length

		// N operation:
		case sam.CigarSkipped:
			exonStart := currBlockStart              // Previous exon starting here.
			exonEnd := currBlockStart + currBlockLen // Previous exon ends here.

			// Create exon object:
			exonId := fmt.Sprintf("exon_%d", exonNr)
			exon := gene.Exon{transcript, exonStart - record.Pos, exonEnd - exonStart, exonId}

			// Discard zero length exons - FIXME: maybe this should not happen.
			if exon.Len() > 0 {
				exons = append(exons, exon) // Register exon.
			}

			currBlockLen = 0                  // Reset exon length counter.
			currBlockStart = exonEnd + length // Next exon starts after the N operation.
			exonNr++

		default:
			L.Fatalf("Unsupported CIGAR operation %s\n in record %s\n", op, record.Name) // FIXME
		}

	}

	// Deal with the last exon:
	exonStart := currBlockStart
	exonEnd := currBlockStart + currBlockLen
	exonId := fmt.Sprintf("exon_%d", exonNr)
	exon := gene.Exon{transcript, exonStart - record.Pos, exonEnd - exonStart, exonId}
	if exon.Len() > 0 {
		exons = append(exons, exon)
	}

	// Add exons to the transcript:
	err := transcript.SetExons(exons...)
	if err != nil {
		L.Fatalf("Could not set exons for %s: %s\n", transcript.ID, err)
	}

	// Convert transcript into GFF2 features:
	trFeatures := Transcript2GFF(transcript)
	return trFeatures
}

// Write a slice of GFF fetures to file.
func WriteFeatures(features []gff.Feature, writer *gff.Writer) {
	// Write GFF features:
	for _, feat := range features {
		_, err := writer.Write(&feat)
		if err != nil {
			L.Fatalf("Failed to write feature %s: %s", feat, err)
		}
	}
}

// Convert a gene.CodingTranscript object into a slice of GFF features.
func Transcript2GFF(tr *gene.CodingTranscript) []gff.Feature {
	res := make([]gff.Feature, 0, len(tr.Exons())+1)

	trFeat := gff.Feature{
		SeqName:        tr.Location().Name(),
		Source:         "pinfish",
		Feature:        "mRNA",
		FeatStart:      tr.Start(),
		FeatEnd:        tr.End(),
		FeatScore:      nil,
		FeatStrand:     seq.Strand(tr.Orient),
		FeatFrame:      gff.NoFrame,
		FeatAttributes: gff.Attributes{gff.Attribute{Tag: "gene_id", Value: "\"" + tr.ID + "\""}, gff.Attribute{Tag: "transcript_id", Value: "\"" + tr.ID + "\";"}},
	}

	res = append(res, trFeat)

	for _, exon := range tr.Exons() {
		exFeat := gff.Feature{
			SeqName:        tr.Location().Name(),
			Source:         "pinfish",
			Feature:        "exon",
			FeatStart:      tr.Offset + exon.Start(),
			FeatEnd:        tr.Offset + exon.End(),
			FeatScore:      nil,
			FeatStrand:     seq.Strand(tr.Orient),
			FeatFrame:      gff.NoFrame,
			FeatAttributes: gff.Attributes{gff.Attribute{Tag: "transcript_id", Value: "\"" + tr.ID + "\";"}},
		}
		res = append(res, exFeat)

	}

	return res
}
