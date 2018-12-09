# Python query engine for PyRice package

## Function

### Query Module

```py
MultiQuery().check_gene(idents=[],locs=[])
```
Return 2 set :id and loc and 2 dict: id-loc and loc-id (in "/support/id_dict" and "/support/loc_dict")
1. idents: list of idents
2. locs: list of locs

```py
MultiQuery().search_gene(chro="", start_pos="",end_pos="",dbs=[])
```
Return a dictionary in snpseek database
1. dbs : list of database,default is search gene in all of "msu7","rap","iric"

```py
multiquery = MultiQuery().
multiquery.query_all(idents=[],locs=[],dbs=[])
```
Return a dictionary in databases follow idents and locs
1. idents: list of idents
2. locs: list of locs
3. dbs : list of database,default is search gene in all of "oryzabase", "Gramene", "funricegene_genekeywords",
                       "funricegene_faminfo", "msu", "rapdb","ic4r",
                       "funricegene_geneinfo"

```py
multiquery = MultiQuery()
multiquery.query_all(idents=[],locs=[],dbs="all")
multiquery.save(folder_path)
```
Save file after query
1. folder_path/data/db.csv: save all gene in all database in one file (.csv)
2. folder_path/geng/...: save each gene in all database (.txt)

## Structure of Database description

```xml
<database dbname="name of the database" type="Type of the response" method="GET or POST">
    <link stern="the link section before the query" aft="section behind the query"/>
    <headers>
        <header type="">Column number 1</header>
        <header type="">Column number 2</header>
        etc.
    </headers>
    <fields>
        <field>Query argument number 1</field>
    </fields>
    <data_struct indicator="indicator of return data segment" identifier="the attribute to identify data section" identification_string="value of said identifier" line_separator="indicator of a line of data" cell_separator="indicator of a cell of data"/>
    <prettify>Regular expression of unwanted character</prettify>
</database>
```

## Example run in Pycharm

### Example of system check gene

```bash
set_ids, set_locs, idents, locs = MultiQuery().check_gene(idents=["Os08g0164400", "Os07g0586200"],locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914"])
```
```bash
Set of ids {'Os07g0586200', 'Os10g0206500', 'Os08g0164400'} 
Set of locs {'LOC_Os08g06740', 'LOC_Os07g39750', 'LOC_Os08g06760'}
```

### Example of system search gene in database snpeek

```bash
 test = MultiQuery().search_gene(chro="chr01", start_pos="1",
                                    end_pos="10000",dbs="all")
```
```bash
Database: snpseek
rap [{'contig': 'chr01', 'fmin': 2982, 'fmax': 10815, 'uniquename': 'Os01g0100100', 'strand': 1, 'msu7Name': 'LOC_Os01g01010', 'raprepName': 'Os01g0100100', 'rappredName': None, 'iricname': 'OsNippo01g010050', 'fgeneshName': 'chr01-gene_1', 'description': 'RabGAP/TBC domain containing protein. (Os01t0100100-01)'}]
iric [{'contig': 'chr01', 'fmin': 2902, 'fmax': 10817, 'uniquename': 'OsNippo01g010050', 'strand': 1, 'msu7Name': 'LOC_Os01g01010', 'raprepName': 'Os01g0100100', 'rappredName': None, 'iricname': 'OsNippo01g010050', 'fgeneshName': 'chr01-gene_1', 'description': 'RabGAP/TBC domain containing protein. (Os01t0100100-01)'}]
msu7 [{'contig': 'chr01', 'fmin': 2902, 'fmax': 10817, 'uniquename': 'LOC_Os01g01010', 'strand': 1, 'msu7Name': 'LOC_Os01g01010', 'raprepName': 'Os01g0100100', 'rappredName': None, 'iricname': 'OsNippo01g010050', 'fgeneshName': 'chr01-gene_1', 'description': 'TBC domain containing protein, expressed'}]
```

### Example of query all database 

```bash
multi_query = MultiQuery()
test = multi_query.query_all(idents=["Os08g0164400", "Os07g0586200","Os01g0100900"],locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914","LOC_Os01g01019"],dbs='all')
#save_file
multi_query.save_file("./result/")
```
One gene (more detail in file)
```bash
Gene: Os01g0100200-LOC_Os01g01019
Gramene [{'_id': 'Os01g0100200', 'name': 'Os01g0100200', 'description': 'Conserved hypothetical protein. (Os01t0100200-01)', 'biotype': 'protein_coding', 'taxon_id': 39947, 'system_name': 'oryza_sativa', 'db_type': 'core', 'gene_idx': 2, 'location': {'region': '1', 'start': 11218, 'end': 12435, 'strand': 1, 'map': 'GCA_001433935.1'}, 'xrefs': [{'db': 'UniParc', 'ids': ['UPI000043A2EB']}, {'db': 'Uniprot/SPTREMBL', 'ids': ['Q655L9']}, {'db': 'protein_id', 'ids': ['EEE53690.1', 'BAD45493.1', 'BAH90846.1', 'BAE79747.1', 'BAS69910.1']}], 'gene_structure': {'exons': [{'id': 'Os01t0100200-01.exon1', 'start': 1, 'end': 843}, {'id': 'Os01t0100200-01.exon2', 'start': 935, 'end': 1218}], 'transcripts': [{'exons': ['Os01t0100200-01.exon1', 'Os01t0100200-01.exon2'], 'length': 1127, 'exon_junctions': [843], 'cds': {'start': 581, 'end': 1009}, 'translation': {'id': 'Os01t0100200-01', 'length': 142, 'features': {}}, 'id': 'Os01t0100200-01'}], 'canonical_transcript': 'Os01t0100200-01'}, 'annotations': {'taxonomy': {'entries': [{'_id': 39947, 'name': 'Oryza sativa Japonica Group'}], 'ancestors': [1, 2759, 3193, 3398, 4447, 4479, 4527, 4530, 4734, 33090, 35493, 38820, 58023, 58024, 78536, 131221, 131567, 147367, 147380, 359160, 1437183, 1437197, 1648021]}, 'familyRoot': {'entries': [{'_id': 4527, 'name': 'Oryza'}], 'ancestors': [1, 2759, 3193, 3398, 4447, 4479, 4734, 33090, 35493, 38820, 58023, 58024, 78536, 131221, 131567, 147367, 147380, 359160, 1437183, 1437197, 1648021]}}, 'homology': {'gene_tree': {'id': 'EPlGT00140000023342', 'root_taxon_id': 4527, 'root_taxon_name': 'Oryza', 'duplications': [40149]}, 'homologous_genes': {'ortholog_one2one': ['BGIOSGA002570', 'ONIVA02G14150', 'OBART01G00030', 'OB02G44720'], 'ortholog_one2many': ['OMERI05G17330', 'OMERI06G01950'], 'syntenic_ortholog_one2one': ['OGLUM01G00040']}}, 'bins': {'fixed_100': 3106, 'fixed_200': 6206, 'fixed_500': 15506, 'fixed_1000': 31008, 'uniform_1Mb': 31350, 'uniform_2Mb': 15789, 'uniform_5Mb': 6466, 'uniform_10Mb': 3342}, 'species_idx': 5}]
rapdb {'ID': 'Os01g0100200', 'Description': 'Conserved hypothetical protein. (Os01t0100200-01)', 'Position': 'chr01:11218..12435', 'RAP-DB Gene Symbol Synonym(s)': '', 'RAP-DB Gene Name Synonym(s)': '', 'CGSNL Gene Symbol': '', 'CGSNL Gene Name': '', 'Oryzabase Gene Symbol Synonym(s)': '', 'Oryzabase Gene Name Synonym(s)': ''}
ic4r {'All': '', 'OS Gene ID': 'Os01g0100200', 'LOC Gene ID': 'LOC_Os01g01019', 'Symbol': '-', 'Location': 'Chr1:11218-12435 (+)', 'Description': 'expressed
```

### Example of database description

```xml
<database dbname="oryzabase" type="text/html" method="POST">
    <link stern="https://shigen.nig.ac.jp/rice/oryzabase/gene/advanced/list"/>
    <headers>
        <header type="">CGSNL Gene Symbol</header>
        <header type="">Gene symbol synonym(s)</header>
        <header type="">CGSNL Gene Name</header>
        <header type="">Gene name synonym(s)</header>
        <header type="">Chr. No.</header>
        <header type="">Trait Class</header>
        <header type="">Gene Ontology</header>
        <header type="">Trait Ontology</header>
        <header type="">Plant Ontology</header>
        <header type="">RAP ID</header>
        <header type="">Mutant Image</header>
    </headers>
    <fields>
        <field>rapId</field>
    </fields>
    <data_struct indicator="table" identifier="class" identification_string="table_summery_list table_nowrapTh max_width_element" line_separator="tr" cell_separator="td"/>
    <prettify>\n>LOC_.*\n|\n|\r|\t</prettify>
</database>
```

## List of supported database

* Oryzabase
* RapDB
* Gramene
* ic4r
* plntfdb
* SNP-Seek
* funricegene
* MSU
* RiceNetDb
* Uniprot - Get protein ID
* pfam - offline
* Kegg

## List of exception

* Server Exception

    Throw when server response code is not 200.

    Throw with the corresponding server response code.
* Internet Connection Exceptioin

    Throw requests.exceptions.RequestException

    *requests* module exception.
* Timeout Exception

    Throw requests.exceptions.Timeout

    *requests* module exception.
* Database Exception

    Throw when database description is not found.
