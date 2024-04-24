#!/usr/bin/python3
import pprint
from BCBio import GFF
from collections.abc import Generator
from typing import Callable, Self


class GFFRecord:

    def __init__(self,
                 seqid: str = None,
                 source: str = None,
                 type: str = None,
                 start: int = None,
                 end: int = None,
                 score: float = None,
                 strand: chr = None,
                 phase: str = None,
                 attributes: dict[str, str] = None):
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

    @staticmethod
    def parse_line(line: str) -> 'GFFRecord':
        line_split = line.strip().split(sep='\t')

        return GFFRecord(
            seqid=line_split[0] if line_split[0] != "." else None,
            source=line_split[1] if line_split[1] != "." else None,
            type=line_split[2] if line_split[2] != "." else None,
            start=int(line_split[3]) if line_split[3] != "." else None,
            end=int(line_split[4]) if line_split[4] != "." else None,
            score=float(line_split[5]) if line_split[5] != "." else None,
            strand=line_split[6] if line_split[6] != "." else None,
            phase=line_split[7] if line_split[7] != "." else None,
            attributes=GFFRecord._parse_attributes(line_split[8]) if line_split[8] != "." else None,
        )

    @staticmethod
    def _parse_attributes(attributes: str) -> dict[str: str]:
        attributes_dict = dict()
        for attribute in attributes.split(';'):
            key, value = attribute.split('=')
            attributes_dict[key.lower()] = value
        return attributes_dict

    def get_attribute(self, attribute_key: str) -> str:
        return self.attributes.get(attribute_key.lower())

    def get_attribute_dict_value(self, attribute_key: str, attribute_parsed_value_key: str) -> str:
        parsed_value_dict: dict[str, str] = dict()
        attribute_value = self.attributes.get(attribute_key.lower())
        if attribute_value is not None:
            for attribute in attribute_value.split(','):
                key, value = attribute.split(':')
                parsed_value_dict[key.lower()] = value
            return parsed_value_dict.get(attribute_parsed_value_key.lower())

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f"seqid: {self.seqid}\n"
            f"type: {self.type}, source: {self.source}\n"
            f"location: {self.start} - {self.end}\n"
            f"strand: {self.strand}, phase: {self.phase}\n"
            f"attributes: {self.attributes}\n"
            "" if self.score is None else f"score: {self.score}\n"

        )


class GFFParser:

    def __init__(self, generator: Generator[GFFRecord]):
        self.generator = generator

    @staticmethod
    def parse_file(gff_file: str) -> 'GFFParser':
        return GFFParser(GFFParser.__parse_gff3_iterator(gff_file))

    # @staticmethod
    # def combine_parser(parsers: list['GFFParser']) -> 'GFFParser':
    #     return GFFParser(GFFParser.__combine_parser_iterator(parsers))

    # @staticmethod
    # def __combine_parser_iterator(parsers: list['GFFParser']) -> Generator[GFFRecord]:
    #     for parser in parsers:
    #         yield from parser.generator

    @staticmethod
    def __parse_gff3_iterator(file_path: str) -> Generator[GFFRecord]:
        # writing my own parser
        with open(file_path) as in_handle:
            for line in in_handle:
                if line.startswith('#'):
                    continue
                else:
                    yield GFFRecord.parse_line(line)

    def filter(self,
               seqid: str = None,
               source: str = None,
               type: [str] = None,
               start: int = None,
               end: int = None,
               strand: chr = None
               ) -> Self:
        return GFFParser(self.__filter(seqid, source, type, start, end, strand))

    def __filter(self,
                 seqid: str = None,
                 source: str = None,
                 type: [str] = None,
                 start: int = None,
                 end: int = None,
                 strand: chr = None) -> Generator[GFFRecord]:
        for record in self.generator:
            if seqid is not None and record.seqid != seqid: continue
            elif source is not None and record.source != source: continue
            elif type is not None and record.type not in type: continue
            elif start is not None and record.start < start: continue
            # not included < start <= included
                # if record.strand == "+":
                #     if record.start < start: continue
                # elif record.strand == "-":
                #     if record.start > start: continue
            elif end is not None and record.end > end: continue
            # not included <= end < not included
                # if record.strand == "+":
                #     if record.end > end: continue
                # elif record.strand == "-":
                #     if record.end < end: continue
            elif strand is not None and record.strand != strand: continue
            else: yield record

    def filter_gene_id(self, gene_id: int, exclude: bool = False) -> Self:
        return GFFParser(self.__filter_gene_id(gene_id, exclude))

    def __filter_gene_id(self, gene_id: int, exclude: bool = False) -> Generator[GFFRecord]:
        record: GFFRecord
        for record in self.generator:
            record_gene_id = int(record.get_attribute_dict_value('dbxref', 'GeneID'))
            if not exclude and record_gene_id == gene_id: yield record
            elif exclude and record_gene_id != gene_id: yield record
            else: continue

    def filter_attributes(self,
                          attribute_key: str,
                          predicate: Callable[[str], bool]) -> Self:
        return GFFParser(self.__filter_attributes(attribute_key, predicate))

    def __filter_attributes(self,
                            attribute_key: str,
                            predicate: Callable[[str], bool]) -> Generator[GFFRecord]:
        for record in self.generator:
            if (record.attributes.get(attribute_key.lower()) is not None
                    and predicate(record.attributes[attribute_key.lower()])):
                yield record
            else: continue

    def filter_attributes_dict(self,
                               attribute_key: str,
                               attribute_parsed_value_key: str,
                               predicate: Callable[[str], bool]) -> Self:
        def parsed_val_predicate(attribute_value: str) -> bool:
            parsed_value_dict: [str, str] = dict()
            for attribute in attribute_value.split(','):
                key, value = attribute.split(':')
                parsed_value_dict[key.lower()] = value
            return (parsed_value_dict.get(attribute_parsed_value_key.lower()) is not None
                    and predicate(parsed_value_dict[attribute_parsed_value_key.lower()]))

        return self.filter_attributes(attribute_key, parsed_val_predicate)


def examine_gff3(file_path: str):
    examiner = GFF.GFFExaminer()
    with open(file_path) as in_handle:
        pprint.pprint(examiner.available_limits(in_handle))


if __name__ == '__main__':
    examine_gff3("./genomes/7227/ncbi_dataset/data/GCF_000001215.4/genomic.gff")

    for _record in (GFFParser.parse_file("./genomes/7227/ncbi_dataset/data/GCF_000001215.4/genomic.gff")
                    .filter(type=["gene", "exon"])
                    .generator):
        print(_record)

    print("###################")

    for _record in (GFFParser.parse_file("./genomes/7227/ncbi_dataset/data/GCF_000001215.4/genomic.gff")
                    .filter_attributes(attribute_key="parent", predicate=(lambda parent: "gene-Dmel_CG41624" == parent))
                    .generator):
        print(_record)


# GFF3 files are nine-column, tab-delimited, plain text files.
# based on specifications:  https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

# GFF3 files are nine-column, tab-delimited, plain text files. Literal use of tab, newline, carriage return,
# the percent (%) sign, and control characters must be encoded using RFC 3986 Percent-Encoding;
# no other characters may be encoded.
# Backslash and other ad-hoc escaping conventions that have been added to the GFF format are not allowed.
# The file contents may include any character in the set supported by the operating environment,
# although for portability with other systems, use of UTF-8 is recommended.

# encoding:
#   tab (%09)
#   newline (%0A)
#   carriage return (%0D)
#   % percent (%25)
#   control characters (%00 through %1F, %7F)

# Column 1: "seqid"
#     The ID of the landmark used to establish the coordinate system for the current feature.
#     IDs may contain any characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|].
#     In particular, IDs may not contain unescaped whitespace and must not begin with an unescaped ">".
# Column 2: "source"
#     The source is a free text qualifier intended to describe the algorithm
#     or operating procedure that generated this feature.
#     Typically this is the name of a piece of software, such as "Genescan" or a database name, such as "Genbank."
#     In effect, the source is used to extend the feature ontology by adding a qualifier to the
#     type creating a new composite type that is a subclass of the type in the type column.
# Column 3: "type"
#     The type of the feature (previously called the "method").
#     This is constrained to be either a term from the Sequence Ontology or an SO accession number.
#     The latter alternative is distinguished using the syntax SO:000000.
#     In either case, it must be sequence_feature (SO:0000110) or an is_a child of it.
# Columns 4 & 5: "start" and "end"
#     The start and end coordinates of the feature are given in positive 1-based integer coordinates,
#     relative to the landmark given in column one. Start is always less than or equal to end.
#     For features that cross the origin of a circular feature
#     (e.g. most bacterial genomes, plasmids, and some viral genomes),
#     the requirement for start to be less than or equal to end is satisfied by making end = the position
#     of the end + the length of the landmark feature.
#     For zero-length features, such as insertion sites,
#     start equals end and the implied site is to the right of the indicated base in the direction of the landmark.
# Column 6: "score"
#     The score of the feature, a floating point number. As in earlier versions of the format,
#     the semantics of the score are ill-defined.
#     It is strongly recommended that E-values be used for sequence similarity features,
#     and that P-values be used for ab initio gene prediction features.
# Column 7: "strand"
#     The strand of the feature. + for positive strand (relative to the landmark),
#     - for minus strand, and . for features that are not stranded. In addition,
#     ? can be used for features whose strandedness is relevant, but unknown.
# Column 8: "phase"
#     For features of type "CDS", the phase indicates where the next codon begins relative to the 5' end
#     (where the 5' end of the CDS is relative to the strand of the CDS feature) of the current CDS feature.
#     For clarification the 5' end for CDS features on the plus strand is the feature's start
#     and and the 5' end for CDS features on the minus strand is the feature's end.
#     The phase is one of the integers 0, 1, or 2, indicating the number of
#     bases forward from the start of the current CDS feature the next codon begins.
#     A phase of "0" indicates that a codon begins on the first nucleotide of the CDS feature (i.e. 0 bases forward),
#     a phase of "1" indicates that the codon begins at the second nucleotide of this CDS feature and a phase of "2"
#     indicates that the codon begins at the third nucleotide of this region.
#     Note that ‘Phase’ in the context of a GFF3 CDS feature should not be confused with the similar concept
#     of frame that is also a common concept in bioinformatics.
#     Frame is generally calculated as a value for a given base relative to the start of the complete open reading
#     frame (ORF) or the codon (e.g. modulo 3) while CDS phase describes
#     the start of the next codon relative to a given CDS feature.
#     The phase is REQUIRED for all CDS features.
# Column 9: "attributes"
#     A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons.
#     URL escaping rules are used for tags or values containing the following characters: ",=;".
#     Spaces are allowed in this field, but tabs must be replaced with the %09 URL escape.
#     Attribute values do not need to be and should not be quoted.
#     The quotes should be included as part of the value by parsers and not stripped.
#
#     These tags have predefined meanings:
#     ID
#         Indicates the ID of the feature. The ID attribute is required for features that have children
#         (e.g. gene and mRNAs), or for those that span multiple lines, but are optional for other features.
#         IDs for each feature must be unique within the scope of the GFF file.
#         In the case of discontinuous features (i.e. a single feature that exists over multiple genomic locations)
#         the same ID may appear on multiple lines.
#         All lines that share an ID must collectively represent a single feature.
#     Name
#         Display name for the feature. This is the name to be displayed to the user. Unlike IDs,
#         there is no requirement that the Name be unique within the file.
#     Alias
#         A secondary name for the feature.
#         It is suggested that this tag be used whenever a secondary identifier for the feature is needed,
#         such as locus names and accession numbers. Unlike ID,
#         there is no requirement that Alias be unique within the file.
#     Parent
#         Indicates the parent of the feature. A parent ID can be used to group exons into transcripts,
#         transcripts into genes, an so forth. A feature may have multiple parents.
#         Parent can only be used to indicate a partof relationship.
#     Target
#         Indicates the target of a nucleotide-to-nucleotide or protein-to-nucleotide alignment.
#         The format of the value is "target_id start end [strand]", where strand is optional and may be "+" or "-".
#         If the target_id contains spaces, they must be escaped as hex escape %20.
#     Gap
#         The alignment of the feature to the target if the two are not collinear (e.g. contain gaps).
#         The alignment format is inspired from the CIGAR format described in the Exonerate documentation.
#     Derives_from
#         Used to disambiguate the relationship between one feature and another when the relationship is a temporal
#         one rather than a purely structural "part of" one.
#         This is needed for polycistronic genes.
#         See "PATHOLOGICAL CASES" for further discussion.
#     Note
#         A free text note.
#     Dbxref
#         A database cross reference.
#         See the section "Ontology Associations and Db Cross References" for details on the format.
#     Ontology_term
#         A cross reference to an ontology term.
#         See the section "Ontology Associations and Db Cross References" for details.
#     Is_circular
#         A flag to indicate whether a feature is circular.
#         See extended discussion below.
#
#     Multiple attributes of the same type are indicated by separating the values with the comma "," character, as in:
#     Parent=AF2312,AB2812,abc-3
#
#     In addition to Parent, the Alias, Note, Dbxref and Ontology_term attributes can have multiple values.
#     Note that attribute names are case sensitive. "Parent" is not the same as "parent".
#
#     All attributes that begin with an uppercase letter are reserved for later use.
#     Attributes that begin with a lowercase letter can be used freely by applications.

# SeqRecord:
#    seqid – SeqRecord ID
#    source – Feature qualifier with key “source”
#    type – Feature type attribute
#    start, end – The Feature Location
#    score – Feature qualifier with key “score”
#    strand – Feature strand attribute
#    phase – Feature qualifier with key “phase”

# sample
# NT_033777.3	RefSeq	region	1	32079331	.	+	.	ID=NT_033777.3:1..32079331;Dbxref=taxon:7227;Name=3R;chromosome=3R;gbkey=Src;genome=chromosome;genotype=y[1]%3B Gr22b[1] Gr22d[1] cn[1] CG33964[R4.2] bw[1] sp[1]%3B LysC[1] MstProx[1] GstD5[1] Rh6[1];mol_type=genomic DNA
# NT_033777.3	RefSeq	mobile_genetic_element	97163	98229	.	-	.	ID=id-NT_033777.3:97163..98229;Dbxref=FLYBASE:FBti0215171;gbkey=mobile_element;mobile_element_type=transposon:1360{}6392
# NT_033777.3	RefSeq	mobile_genetic_element	379801	384688	.	+	.	ID=id-NT_033777.3:379801..384688;Dbxref=FLYBASE:FBti0215319;gbkey=mobile_element;mobile_element_type=transposon:invader2{}6541
# NT_033777.3	RefSeq	mobile_genetic_element	410857	411958	.	+	.	ID=id-NT_033777.3:410857..411958;Dbxref=FLYBASE:FBti0215172;gbkey=mobile_element;mobile_element_type=transposon:1360{}6393
# NT_033777.3	RefSeq	mobile_genetic_element	418655	419755	.	+	.	ID=id-NT_033777.3:418655..419755;Dbxref=FLYBASE:FBti0215173;gbkey=mobile_element;mobile_element_type=transposon:1360{}6394
# NT_033777.3	RefSeq	mobile_genetic_element	420748	421848	.	+	.	ID=id-NT_033777.3:420748..421848;Dbxref=FLYBASE:FBti0215174;gbkey=mobile_element;mobile_element_type=transposon:1360{}6395
# NT_033777.3	RefSeq	mobile_genetic_element	533820	534896	.	+	.	ID=id-NT_033777.3:533820..534896;Dbxref=FLYBASE:FBti0215175;gbkey=mobile_element;mobile_element_type=transposon:1360{}6396
# NT_033777.3	RefSeq	mobile_genetic_element	562645	563742	.	+	.	ID=id-NT_033777.3:562645..563742;Dbxref=FLYBASE:FBti0215176;gbkey=mobile_element;mobile_element_type=transposon:1360{}6397
# NT_033777.3	RefSeq	gene	567076	2532932	.	+	.	ID=gene-Dmel_CG45784;Dbxref=FLYBASE:FBgn0267431,GeneID:26067052;Name=Myo81F;Note=Annotated by Drosophila Heterochromatin Genome Project%2C Lawrence Berkeley National Lab%2C http://www.dhgp.org;cyt_map=81F1-81F3;description=Myosin 81F;gbkey=Gene;gen_map=3-47.1 cM;gene=Myo81F;gene_biotype=protein_coding;gene_synonym=CG15831,CG40155,CG40204,CG41281,CG41518,CG41527,CG42621,CG42622,CG42623,CG42624,CG45784,Dmel\CG45784;locus_tag=Dmel_CG45784;old_locus_tag=Dmel_CG15831%2CDmel_CG40155%2CDmel_CG40204%2CDmel_CG41281%2CDmel_CG41518%2CDmel_CG41527%2CDmel_CG42621%2CDmel_CG42622%2CDmel_CG42623%2CDmel_CG42624
# NT_033777.3	RefSeq	mRNA	567076	2532932	.	+	.	ID=rna-NM_001316563.1;Parent=gene-Dmel_CG45784;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Name=NM_001316563.1;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	567076	567268	.	+	.	ID=exon-NM_001316563.1-1;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	835376	835491	.	+	.	ID=exon-NM_001316563.1-2;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	869486	869548	.	+	.	ID=exon-NM_001316563.1-3;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
# NT_033777.3	RefSeq	exon	895786	895893	.	+	.	ID=exon-NM_001316563.1-4;Parent=rna-NM_001316563.1;Dbxref=FLYBASE:FBtr0392909,GeneID:26067052,GenBank:NM_001316563.1,FLYBASE:FBgn0267431;Note=Myo81F-RB%3B Dmel\Myo81F-RB%3B CG45784-RB%3B Dmel\CG45784-RB;gbkey=mRNA;gene=Myo81F;locus_tag=Dmel_CG45784;orig_protein_id=gnl|FlyBase|CG45784-PB|gb|ALI30523;orig_transcript_id=gnl|FlyBase|CG45784-RB;product=myosin 81F;transcript_id=NM_001316563.1
