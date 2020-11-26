package main

import (
	"bufio"
	"github.com/biogo/hts/bam"
	"os"
)

// Create new BAM reader from file.
func NewBamReader(bamFile string, nrProc int) *bam.Reader {
	fh := os.Stdin
	var err error
	if bamFile != "_" {
		fh, err = os.Open(bamFile)
		if err != nil {
			L.Fatalf("Could not open input file %s: %s\n", bamFile, err)
		}
	}

	reader, err := bam.NewReader(bufio.NewReader(fh), nrProc)
	if err != nil {
		L.Fatalf("Could not create BAM reader for %s: %s\n", bamFile, err)
	}
	return reader
}
