import pysam

reverse_strand = {0: 16, 16: 0}

with pysam.AlignmentFile(snakemake.input['bam'], "rb") as input_bam:
    with pysam.AlignmentFile(snakemake.output['sam_reverse'], "w", template=input_bam) as output_sam:

        for read in input_bam:
            if read.has_tag('ts') and read.flag in reverse_strand:
                if read.get_tag('ts') == '-':
                    read.flag = reverse_strand[read.flag]
                    read.set_tag('ts', '+')

            output_sam.write(read)
