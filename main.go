package main

import (
	"log"
	"os"
	"runtime"
)

func main() {
	L = NewLogger("spliced_bam2gff: ", log.Ltime)

	// Parse command line arguments:
	args := new(CmdArgs)
	args.Parse()

	// Set the maximum number of OS threads to use:
	runtime.GOMAXPROCS(int(args.MaxProcs))

	if args.OutDir == "" {
		// Convert spliced BAM entries to GFF transcripts:
		SplicedBam2GFF(args.InputFiles[0], os.Stdout, int(args.MaxProcs), args.MinimapInput, args.StrandBehaviour, int(args.MaxDel), args.KeepS)
	} else {
		SplicedBam2PartGFF(args.InputFiles[0], args.OutDir, int(args.MinBundle), int(args.MaxProcs), args.MinimapInput, args.StrandBehaviour, int(args.MaxDel), args.KeepS)

	}

}
