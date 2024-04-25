#!/usr/bin/python3
import re
import sys

from gff_parser import GFFRecord, GFFParser
from collections import defaultdict
import Bio.SeqIO.FastaIO as FastaIO
# from Bio.Seq import Seq


class SequenceRecord:

    # (gff_record, sequence)
    def __init__(self, parent: 'SequenceRecordsBundle', gff_record: GFFRecord, sequence: str):
        self.parent = parent
        self.gff_record = gff_record
        self.sequence = sequence

    def get_taxonomy_id(self) -> int: return self.parent.taxonomy_id

    def get_genome_id(self) -> str: return self.parent.genome_id

    @staticmethod
    def from_genome_sequence(parent: 'SequenceRecordsBundle',
                             gff_record: GFFRecord,
                             genome_sequence: str) -> 'SequenceRecord':
        # counting starts from 1
        # get slice from genome sequence according to start/end in gff3 record

        gene_sequence = genome_sequence[gff_record.start-1:gff_record.end-1]
        # if gff_record.strand == "-": gene_sequence = str(Seq(gene_sequence).reverse_complement())

        return SequenceRecord(parent=parent, gff_record=gff_record, sequence=gene_sequence)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f"gff_record: {self.gff_record}"
            f"sequence: {self.sequence}"
        )

    # __gff_attributes_pattern = re.compile(r"gff_attribute:\w+(?=[,; \n\t])")
    __gff_attributes_pattern = re.compile(r"gff_attribute:\w+")

    def as_table_record(self, data_format: str) -> str:

        data_format = data_format.lower()  # yeah... there might be a more elegant solution to this

        # sequence record bundle
        if "sr_taxonomy_id" in data_format:
            data_format = data_format.replace("sr_taxonomy_id", str(self.parent.taxonomy_id))
        if "sr_genome_id" in data_format:
            data_format = data_format.replace("sr_genome_id", str(self.parent.genome_id))
        if "sr_headline" in data_format:
            data_format = data_format.replace("sr_headline", str(self.parent.headline.replace(',', '')))
        if "sr_genome_sequence" in data_format:
            data_format = data_format.replace("sr_genome_sequence", str(self.parent.genome_sequence))

        # sequence record
        if "sr_sequence" in data_format:
            data_format = data_format.replace("sr_sequence", str(self.sequence))

        # gff record
        if "gff_seqid" in data_format:
            data_format = data_format.replace("gff_seqid", str(self.gff_record.seqid))
        if "gff_source" in data_format:
            data_format = data_format.replace("gff_source", str(self.gff_record.source))
        if "gff_type" in data_format:
            data_format = data_format.replace("gff_type", str(self.gff_record.type))
        if "gff_start" in data_format:
            data_format = data_format.replace("gff_start", str(self.gff_record.start))
        if "gff_end" in data_format:
            data_format = data_format.replace("gff_end", str(self.gff_record.end))
        if "gff_score" in data_format:
            data_format = data_format.replace("gff_score", str(self.gff_record.score))
        if "gff_strand" in data_format:
            data_format = data_format.replace("gff_strand", str(self.gff_record.strand))
        if "gff_phase" in data_format:
            data_format = data_format.replace("gff_phase", str(self.gff_record.phase))

        if "gff_gene_id" in data_format:
            data_format = data_format.replace("gff_gene_id", str(
                self.gff_record.get_attribute_dict_value('dbxref', 'GeneID')))

        # if "gff_gene_name" in data_format:
        #     data_format = data_format.replace("gff_gene_name", str(
        #         self.gff_record.get_attribute_dict_value('dbxref', 'Name')))

        shift = 0
        for match in SequenceRecord.__gff_attributes_pattern.finditer(data_format):
            # print(f"match: {match}")
            attribute = match.group(0).split(':', maxsplit=1)[1]
            # print(f"attribute: {attribute}")
            attribute_value = self.gff_record.get_attribute(attribute)
            # print(f"attribute_value: {attribute_value}")
            # print(f"attributes: {self.gff_record.attributes}")
            if attribute_value is None: sys.stderr.write(f"No value found for attribute {attribute}\n")
            # print(f"shift: {shift}")
            # print(f"line:\n\t{data_format[:match.start(0)+shift]} + {str(attribute_value)} + {data_format[match.end(0)+shift:]}")
            # print()
            data_format = data_format[:match.start(0)+shift] + str(attribute_value) + data_format[match.end(0)+shift:]
            shift += len(str(attribute_value)) - (match.end(0) - match.start(0))

        return data_format


class SequenceRecordsBundle:

    # ((headline, genome_sequence), [(gff_record, sequence)])
    def __init__(self,
                 taxonomy_id: int,
                 genome_id: str,
                 headline: str,
                 genome_sequence: str,
                 sequence_records: list[SequenceRecord]):
        self.taxonomy_id = taxonomy_id
        self.genome_id = genome_id
        self.headline = headline
        self.genome_sequence = genome_sequence  # one genome sequence from seqid as reference
        self.sequence_records = sequence_records

    @staticmethod
    def of_genes(taxonomy_id: int,
                 genome_id: str,
                 gff_file: str,
                 genome_fasta_file: str,
                 filter_types: list[str],
                 gene_id: int) -> list['SequenceRecordsBundle']:
        parser = SequenceRecordsBundle.__create_gff_parser_for_gene_id(gff_file, filter_types, gene_id)
        seqid_record_dict = SequenceRecordsBundle.__collect_to_seqid_record_dict(parser)
        return SequenceRecordsBundle.__parse_from_fasta(taxonomy_id, genome_id, genome_fasta_file, seqid_record_dict)

    @staticmethod
    def of_gene_neighbours(taxonomy_id: int,
                 genome_id: str,
                 gff_file: str,
                 genome_fasta_file: str,
                 filter_types: list[str],
                 gene_id: int,
                 region_width: int) -> list['SequenceRecordsBundle']:
        parser = SequenceRecordsBundle.__create_gff_parser_for_gene_id(gff_file, filter_types, gene_id)
        # TODO rewrite this simpler!
        # parser_combined = GFFParser.combine_parser([SequenceRecordsBundle.__create_gff_parser_for_region(
        #     gff_file, filter_types, (gff_r.start-region_width, gff_r.end+region_width))
        #     for gff_r in parser.generator])
        # seqid_record_dict = SequenceRecordsBundle.__collect_to_seqid_record_dict(parser_combined)
        _sequence_records_bundles_list: list['SequenceRecordsBundle'] = list()
        gff_record: GFFRecord
        for gff_record in parser.generator:
            region_parser = (SequenceRecordsBundle
                             .__create_gff_parser_for_region(gff_file=gff_file, filter_types=filter_types,
                                                             gff_record=gff_record, region_width=region_width))
            seqid_record_dict = SequenceRecordsBundle.__collect_to_seqid_record_dict(region_parser)
            _sequence_records_bundles_list.extend(
                SequenceRecordsBundle.__parse_from_fasta(taxonomy_id, genome_id, genome_fasta_file, seqid_record_dict))
        return _sequence_records_bundles_list

    @staticmethod
    def __create_gff_parser_for_gene_id(gff_file: str,
                                        filter_types: list[str],
                                        gene_id: int) -> GFFParser:
        # gene_ids = [str(gene_id) for gene_id in gene_ids]
        gene_id = str(gene_id)
        return (GFFParser.parse_file(gff_file)  # getting the gff3 records that I need
                .filter(type=filter_types)
                .filter_attributes_dict(attribute_key="Dbxref",
                                        attribute_parsed_value_key="GeneID",
                                        predicate=lambda attribute_parsed_val: attribute_parsed_val == gene_id))
                                            # any([attribute_parsed_val == gene_id for gene_id in gene_ids]))
    @staticmethod
    def __create_gff_parser_for_region(gff_file: str,
                                       filter_types: list[str],
                                       gff_record: GFFRecord,
                                       region_width: int) -> GFFParser:
        # gene_ids = [str(gene_id) for gene_id in gene_ids]
        return (GFFParser.parse_file(gff_file)  # getting the gff3 records that I need
                .filter(type=filter_types,
                        seqid=gff_record.seqid,
                        start=gff_record.start - region_width,
                        end=gff_record.end + region_width)
                .filter_gene_id(gene_id=int(gff_record.get_attribute_dict_value('dbxref', 'GeneID')), exclude=True)
                )
        # any([attribute_parsed_val == gene_id for gene_id in gene_ids])))

    @staticmethod
    def __collect_to_seqid_record_dict(parser: GFFParser) -> dict[str, list[GFFRecord]]:
        seqid_record_dict: dict[str, list[GFFRecord]] = defaultdict(list)  # seqid -> list[records]
        for record in parser.generator:
            seqid_record_dict[record.seqid].append(record)
        return seqid_record_dict

    @staticmethod
    def __parse_from_fasta(taxonomy_id: int, genome_id: str, genome_fasta_file: str,
                           seqid_record_dict: dict[str, list[GFFRecord]]) -> list['SequenceRecordsBundle']:
        gene_records_full_output: list[SequenceRecordsBundle] = list()
        with open(genome_fasta_file, 'r') as file:
            for headline, genome_sequence in FastaIO.SimpleFastaParser(file):
                for seq_id in seqid_record_dict.keys():  # >NT_033777.3 Drosophila melanogaster chromosome 3R
                    if seq_id in headline:  # filter those that have one of the seqids from records in their headline
                        gene_entries: list[SequenceRecord] = list()
                        sequence_records_bundle = SequenceRecordsBundle(taxonomy_id, genome_id, headline, genome_sequence, gene_entries)
                        gene_entries.extend([SequenceRecord.from_genome_sequence(sequence_records_bundle, record, genome_sequence) for record in seqid_record_dict[seq_id]])
                        gene_records_full_output.append(sequence_records_bundle)
        return gene_records_full_output

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f"headline: {self.headline}\n"
            # f"genome_sequence: {self.genome_sequence}\n"
            f"gene_records: {self.sequence_records}\n")


if __name__ == '__main__':

    for _sequence_records_bundle in SequenceRecordsBundle.of_genes(
            taxonomy_id=7227,
            genome_id="GCF_000001215.4",
            gff_file="./genomes/7227/ncbi_dataset/data/GCF_000001215.4/genomic.gff",
            genome_fasta_file="./genomes/7227/ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna",
            filter_types=["gene", "exon"],
            gene_id=41615):
        for _gene_record in _sequence_records_bundle.sequence_records:
            print(f"record: {_gene_record.gff_record}")
            print(f"sequence: {_gene_record.sequence}")


# seq_id is for the specific genomic region
# ##sequence-region NT_033777.3 1 32079331
# NT_033777.3	RefSeq	region	1	32079331	.	+	.	ID=NT_033777.3:1..32079331;Dbxref=taxon:7227;Name=3R;chromosome=3R;gbkey=Src;genome=chromosome;genotype=y[1]%3B Gr22b[1] Gr22d[1] cn[1] CG33964[R4.2] bw[1] sp[1]%3B LysC[1] MstProx[1] GstD5[1] Rh6[1];mol_type=genomic DNA
# NT_033777.3	RefSeq	mobile_genetic_element	97163	98229	.	-	.	ID=id-NT_033777.3:97163..98229;Dbxref=FLYBASE:FBti0215171;gbkey=mobile_element;mobile_element_type=transposon:1360{}6392
# NT_033777.3	RefSeq	gene	567076	2532932	.	+	.	ID=gene-Dmel_CG45784;Dbxref=FLYBASE:FBgn0267431,GeneID:26067052;Name=Myo81F;Note=Annotated by Drosophila Heterochromatin Genome Project%2C Lawrence Berkeley National Lab%2C http://www.dhgp.org;cyt_map=81F1-81F3;description=Myosin 81F;gbkey=Gene;gen_map=3-47.1 cM;gene=Myo81F;gene_biotype=protein_coding;gene_synonym=CG15831,CG40155,CG40204,CG41281,CG41518,CG41527,CG42621,CG42622,CG42623,CG42624,CG45784,Dmel\CG45784;locus_tag=Dmel_CG45784;old_locus_tag=Dmel_CG15831%2CDmel_CG40155%2CDmel_CG40204%2CDmel_CG41281%2CDmel_CG41518%2CDmel_CG41527%2CDmel_CG42621%2CDmel_CG42622%2CDmel_CG42623%2CDmel_CG42624
# NT_033777.3	RefSeq	mRNA	567076	2532932	.	+	.	ID=rna-NM_001316563.1;Parent=gene-Dmel_CG45784;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Name=NM_001316563.1;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	567076	567268	.	+	.	ID=exon-NM_001316563.1-1;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1