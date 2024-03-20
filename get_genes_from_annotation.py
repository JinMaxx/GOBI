#!/usr/bin/python3

from gff3_parser import GFF_Record, GFF3Parser
import Bio.SeqIO.FastaIO as FastaIO


# TODO! It works!! Now do this for the other strand

class GeneRecord:

    # (gff_records, gene_sequence)
    def __init__(self, gff_record: GFF_Record, gene_sequence: str):
        self.gff_record = gff_record
        self.gene_sequence = gene_sequence

    @staticmethod
    def generate_from_genome_sequence(gff_record: GFF_Record, genome_sequence: str) -> 'GeneRecord':
        return GeneRecord(gff_record=gff_record, gene_sequence=genome_sequence[gff_record.start-1:gff_record.end-1])


class GeneRecordFull:

    # ((headline, genome_sequence), [(gff_records, gene_sequence)])
    def __init__(self, headline: str, genome_sequence: str, gene_records: list[GeneRecord]):
        self.headline = headline
        self.genome_sequence = genome_sequence
        self.gene_records = gene_records


def parse_attribute_value_dict(attribute_val: str) -> dict[str, str]:
    # Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431
    # -> {FLYBASE: FBtr0392909, GeneID: 26067052, GenBank: NM_001316563.1, FLYBASE: FBgn0267431}
    attribute_value_dict = dict()
    for attribute in attribute_val.split(','):
        key, value = attribute.split(':')
        attribute_value_dict[key] = value
    # if attribute_value_dict.get("GeneID") == str(41615): print(attribute_value_dict)  # TEST
    return attribute_value_dict


def _create_parser_for_gene_id(gff_file: str, gene_id: int, filter_types: list[str]) -> GFF3Parser:
    attribute_value_includes_gene_id_predicate = lambda attribute_value: parse_attribute_value_dict(attribute_value).get("GeneID") == str(gene_id)

    # getting the gff3 records that I need
    return (GFF3Parser.parse_file(gff_file)
            .filter(type=filter_types)
            .filter_attributes(attribute_key="Dbxref", predicate=attribute_value_includes_gene_id_predicate))


def _transform_to_seq_id_to_record_dict(parser: GFF3Parser) -> dict[str, list[GFF_Record]]:
    # NT_033779.5
    # ->
    # seqid: NT_033779.5
    # source RefSeq
    # type exon
    # start 18538030
    # end 18538185
    # score None
    # strand: -
    # phase: None
    # attributes: {'id': 'exon-NM_136041.2-8', 'parent': 'rna-NM_136041.2',
    #   'dbxref': 'FLYBASE:FBtr0112544,GeneID:35107,GenBank:NM_136041.2,FLYBASE:FBgn0085370',
    #   'note': 'Pde11-RB%3B Dmel\\Pde11-RB%3B CG34341-RB%3B Dmel\\CG34341-RB',
    #   'gbkey': 'mRNA', 'gene': 'Pde11', 'locus_tag': 'Dmel_CG34341',
    #   'orig_protein_id': 'gnl|FlyBase|CG34341-PB|gb|AAF53675', 'orig_transcript_id': 'gnl|FlyBase|CG34341-RB',
    #   'product': 'phosphodiesterase 11%2C transcript variant B', 'transcript_id': 'NM_136041.2'}

    seq_id_to_record_dict: dict[str, list[GFF_Record]] = dict()  # seqid -> list[records]
    record: GFF_Record  # type annotation
    for record in parser.generator:
        if seq_id_to_record_dict.get(record.seqid) is None:
            seq_id_to_record_dict[record.seqid] = []  # create new list
        seq_id_to_record_dict[record.seqid].append(record)
    return seq_id_to_record_dict

def _generate_gene_records_full(genome_fasta_file: str,
                                seq_id_to_record_dict: dict[str, list[GFF_Record]]) -> list[GeneRecordFull]:
    gene_records_full: [GeneRecordFull] = list()  # [((headline, genome_sequence), [(gff_records, gene_sequence)])]
    # while iterating over all sequences get those that have one of those seqids from the records in their headline

    with open(genome_fasta_file, 'r') as file:
        for headline, genome_sequence in FastaIO.SimpleFastaParser(file):
            # if any(seq_id in headline for seq_id in seq_id_to_record_dict.keys()):
            for seq_id in seq_id_to_record_dict.keys():
                if seq_id in headline:  # >NT_033777.3 Drosophila melanogaster chromosome 3R
                    # get the slice of gene from genome sequence according to seqid
                    # counting starts from 1
                    gene_entries: list[GeneRecord] = [GeneRecord.generate_from_genome_sequence(record, genome_sequence)
                                                      for record in seq_id_to_record_dict[seq_id]]
                    # append them with one genome sequence as reference
                    gene_records_full.append(GeneRecordFull(headline, genome_sequence, gene_entries))

    return gene_records_full


# TODO: maybe write a class for that
def get_records_for_genes(gff_file: str, genome_fasta_file: str, gene_id: int) -> [GeneRecordFull]:
    parser = _create_parser_for_gene_id(gff_file, gene_id, filter_types=["gene", "exon"])
    seq_id_to_record_dict = _transform_to_seq_id_to_record_dict(parser)
    gene_records_full = _generate_gene_records_full(genome_fasta_file, seq_id_to_record_dict)

    return gene_records_full


if __name__ == '__main__':
    gff_file = "./genomes/7227/ncbi_dataset/data/GCF_000001215.4/genomic.gff"
    fasta_file = "./genomes/7227/ncbi_dataset/data/GCF_000001215.4/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
    gene_id = 41615

    for gene_record_full in get_records_for_genes(gff_file, fasta_file, gene_id):
        for gene_record in gene_record_full.gene_records:
            print(f"record: {gene_record.gff_record}")
            print(f"sequence: {gene_record.gene_sequence}")


# seq_id is for the specific genomic region

# ##sequence-region NT_033777.3 1 32079331
# ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7227
# NT_033777.3	RefSeq	region	1	32079331	.	+	.	ID=NT_033777.3:1..32079331;Dbxref=taxon:7227;Name=3R;chromosome=3R;gbkey=Src;genome=chromosome;genotype=y[1]%3B Gr22b[1] Gr22d[1] cn[1] CG33964[R4.2] bw[1] sp[1]%3B LysC[1] MstProx[1] GstD5[1] Rh6[1];mol_type=genomic DNA
# NT_033777.3	RefSeq	mobile_genetic_element	97163	98229	.	-	.	ID=id-NT_033777.3:97163..98229;Dbxref=FLYBASE:FBti0215171;gbkey=mobile_element;mobile_element_type=transposon:1360{}6392
# NT_033777.3	RefSeq	gene	567076	2532932	.	+	.	ID=gene-Dmel_CG45784;Dbxref=FLYBASE:FBgn0267431,GeneID:26067052;Name=Myo81F;Note=Annotated by Drosophila Heterochromatin Genome Project%2C Lawrence Berkeley National Lab%2C http://www.dhgp.org;cyt_map=81F1-81F3;description=Myosin 81F;gbkey=Gene;gen_map=3-47.1 cM;gene=Myo81F;gene_biotype=protein_coding;gene_synonym=CG15831,CG40155,CG40204,CG41281,CG41518,CG41527,CG42621,CG42622,CG42623,CG42624,CG45784,Dmel\CG45784;locus_tag=Dmel_CG45784;old_locus_tag=Dmel_CG15831%2CDmel_CG40155%2CDmel_CG40204%2CDmel_CG41281%2CDmel_CG41518%2CDmel_CG41527%2CDmel_CG42621%2CDmel_CG42622%2CDmel_CG42623%2CDmel_CG42624
# NT_033777.3	RefSeq	mRNA	567076	2532932	.	+	.	ID=rna-NM_001316563.1;Parent=gene-Dmel_CG45784;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Name=NM_001316563.1;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	567076	567268	.	+	.	ID=exon-NM_001316563.1-1;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
