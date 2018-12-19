# Python query engine for PyRice package

## Function

### Query Module


```py
MultiQuery().search_gene(chro="", start_pos="",end_pos="",dbs='all',save_path=None)
```
Return a dictionary in snpseek database
1. dbs : list of database,default is search gene in all: "msu7","rap","iric"
2. save_path: link folder save output function, default is None => don't save 

```py
test = MultiQuery()
file_id = test.search_gene(...)
db = test.query_iric(file_id,dbs='all',save_path=None)
```
Return a dictionary in databases follow iricname
1. file_id : result of function search_gene()
2. dbs : list of database,default is search gene in all of "oryzabase", "Gramene", "funricegene_genekeywords",
                       "funricegene_faminfo", "msu", "rapdb","ic4r",
                       "funricegene_geneinfo"
3. save_path: link folder save output function, default is None => don't save 

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

### Example of system search_gene

```bash
t = time.time()
test = MultiQuery()
file_id = test.search_gene(chro="chr01", start_pos="1",
                            end_pos="10000",dbs="all",save_path="./result/")
print("Time for search gene:",time.time() - t)
print("Ouput file",file_id)
```
```bash
Time for search gene: 5.135738134384155
Ouput file {'OsNippo01g010050': {'msu7Name': {'LOC_Os01g01010'}, 'raprepName': {'Os01g0100100'}}}
```

### Example of system query_iric

```bash
#Search gene and query
t = time.time()
test = MultiQuery()
file_id = test.search_gene(chro="chr01", start_pos="1",
                        end_pos="10000",dbs="all",save_path="./result/")
print("Time for search gene:",time.time() - t)
print("Ouput file",file_id)
#Query iric name
t = time.time()
db = test.query_iric(file_id,dbs="all",save_path="./result/")
print("Time for query:", time.time() - t)
print("Output database",db)
```
```bash
Time for search gene: 5.135738134384155
Ouput file {'OsNippo01g010050': {'msu7Name': {'LOC_Os01g01010'}, 'raprepName': {'Os01g0100100'}}}
Query iricname:  OsNippo01g010050
Time for query: 8.683169841766357
Output database {'OsNippo01g010050': {'Gramene': {'_id': 'Os01g0100100', 'name': 'Os01g0100100', 'description': 'RabGAP/TBC domain containing protein. (Os01t0100100-01)', 'biotype': 'protein_coding', 'taxon_id': 39947, 'system_name': 'oryza_sativa', 'db_type': 'core', 'gene_idx': 1, 'location': {'region': '1', 'start': 2983, 'end': 10815, 'strand': 1, 'map': 'GCA_001433935.1'}, 'xrefs': [{'db': 'UniParc', 'ids': ['UPI0001C7F089']}, {'db': 'Uniprot/SPTREMBL', 'ids': ['A0A0P0UX28']}, {'db': 'STRING', 'ids': ['39947.LOC_Os01g01010.1']}, {'db': 'RefSeq_peptide', 'ids': ['XP_015622096.1']}, {'db': 'protein_id', 'ids': ['BAS69908.1']}, {'db': 'RefSeq_dna', 'ids': ['XM_015766610.1']}, {'db': 'EntrezGene', 'ids': ['4326813']}], 'gene_structure': {'exons': [{'id': 'Os01t0100100-01.exon1', 'start': 1, 'end': 286}, {'id': 'Os01t0100100-01.exon2', 'start': 372, 'end': 634}, {'id': 'Os01t0100100-01.exon3', 'start': 1375, 'end': 1473}, {'id': 'Os01t0100100-01.exon4', 'start': 2475, 'end': 2578}, {'id': 'Os01t0100100-01.exon5', 'start': 4154, 'end': 4962}, {'id': 'Os01t0100100-01.exon6', 'start': 5046, 'end': 5168}, {'id': 'Os01t0100100-01.exon7', 'start': 5250, 'end': 5338}, {'id': 'Os01t0100100-01.exon8', 'start': 5426, 'end': 5626}, {'id': 'Os01t0100100-01.exon9', 'start': 6228, 'end': 6633}, {'id': 'Os01t0100100-01.exon10', 'start': 7120, 'end': 7205}, {'id': 'Os01t0100100-01.exon11', 'start': 7292, 'end': 7448}, {'id': 'Os01t0100100-01.exon12', 'start': 7522, 'end': 7833}], 'transcripts': [{'exons': ['Os01t0100100-01.exon1', 'Os01t0100100-01.exon2', 'Os01t0100100-01.exon3', 'Os01t0100100-01.exon4', 'Os01t0100100-01.exon5', 'Os01t0100100-01.exon6', 'Os01t0100100-01.exon7', 'Os01t0100100-01.exon8', 'Os01t0100100-01.exon9', 'Os01t0100100-01.exon10', 'Os01t0100100-01.exon11', 'Os01t0100100-01.exon12'], 'length': 2935, 'exon_junctions': [286, 549, 648, 752, 1561, 1684, 1773, 1974, 2380, 2466, 2623], 'cds': {'start': 382, 'end': 2490}, 'translation': {'id': 'Os01t0100100-01', 'length': 702, 'features': {'domain': {'entries': [{'name': 'PF00566', 'description': 'RabGAP-TBC', 'db': 'Pfam', 'interpro': 'IPR000195', 'start': 389, 'end': 535}, {'name': 'PS50086', 'description': 'TBC_RABGAP', 'db': 'Prosite_profiles', 'interpro': 'IPR000195', 'start': 90, 'end': 537}, {'name': 'SM00164', 'description': 'tbc_4', 'db': 'Smart', 'interpro': 'IPR000195', 'start': 87, 'end': 558}], 'architecture': [{'root': '195', 'start': 87, 'end': 558, 'interpro': 'IPR000195', 'name': 'Rab-GTPase-TBC_dom', 'description': 'Rab-GTPase-TBC domain'}], 'roots': '195'}, 'homologous_superfamily': {'entries': [{'name': 'SSF47923', 'description': 'Ypt/Rab-GAP domain of gyp1p', 'db': 'SuperFamily', 'interpro': 'IPR035969', 'start': 376, 'end': 497}, {'name': 'SSF47923', 'description': 'Ypt/Rab-GAP domain of gyp1p', 'db': 'SuperFamily', 'interpro': 'IPR035969', 'start': 65, 'end': 129}, {'name': 'SSF47923', 'description': 'Ypt/Rab-GAP domain of gyp1p', 'db': 'SuperFamily', 'interpro': 'IPR035969', 'start': 473, 'end': 692}]}}}, 'id': 'Os01t0100100-01'}], 'canonical_transcript': 'Os01t0100100-01'}, 'annotations': {'domains': {'entries': [{'id': 'IPR000195', 'name': 'Rab-GTPase-TBC_dom', 'description': 'Rab-GTPase-TBC domain'}, {'id': 'IPR035969', 'name': 'Rab-GTPase_TBC_sf', 'description': 'Rab-GTPase-TBC domain superfamily'}]}, 'GO': {'entries': [{'id': 'GO:0005096', 'name': 'GTPase activator activity', 'namespace': 'molecular_function', 'def': 'Binds to and increases the activity of a GTPase, an enzyme that catalyzes the hydrolysis of GTP.', 'subset': ['gosubset_prok'], 'evidence_code': 'IBA'}, {'id': 'GO:0017137', 'name': 'Rab GTPase binding', 'namespace': 'molecular_function', 'def': 'Interacting selectively and non-covalently with Rab protein, any member of the Rab subfamily of the Ras superfamily of monomeric GTPases.', 'subset': None, 'evidence_code': 'IBA'}, {'id': 'GO:0006886', 'name': 'intracellular protein transport', 'namespace': 'biological_process', 'def': 'The directed movement of proteins in a cell, including the movement of proteins between specific compartments or structures within a cell, such as organelles of a eukaryotic cell.', 'subset': ['gosubset_prok'], 'evidence_code': 'IBA'}, {'id': 'GO:0031338', 'name': 'regulation of vesicle fusion', 'namespace': 'biological_process', 'def': 'Any process that modulates the frequency, rate or extent of vesicle fusion.', 'subset': None, 'evidence_code': 'IBA'}, {'id': 'GO:0090630', 'name': 'activation of GTPase activity', 'namespace': 'biological_process', 'def': 'Any process that initiates the activity of an inactive GTPase through the replacement of GDP by GTP.', 'subset': None, 'evidence_code': 'IBA'}, {'id': 'GO:0005622', 'name': 'intracellular', 'namespace': 'cellular_component', 'def': 'The living contents of a cell; the matter contained within (but not including) the plasma membrane, usually taken to exclude large vacuoles and masses of secretory or ingested material. In eukaryotes it includes the nucleus and cytoplasm.', 'subset': ['goslim_chembl', 'goslim_generic', 'goslim_metagenomics', 'goslim_plant', 'gosubset_prok'], 'evidence_code': 'IBA'}, {'id': 'GO:0012505', 'name': 'endomembrane system', 'namespace': 'cellular_component', 'def': 'A collection of membranous structures involved in transport within the cell. The main components of the endomembrane system are endoplasmic reticulum, Golgi bodies, vesicles, cell membrane and nuclear envelope. Members of the endomembrane system pass materials through each other or though the use of vesicles.', 'subset': ['goslim_aspergillus', 'goslim_candida', 'goslim_yeast'], 'evidence_code': 'IBA'}], 'ancestors': [3674, 5488, 5515, 5575, 5623, 6810, 8047, 8104, 8150, 15031, 15833, 17016, 19899, 30234, 30695, 31267, 32879, 33036, 33043, 34613, 42886, 43085, 43087, 43547, 44093, 44464, 45184, 46907, 50789, 50790, 50794, 51020, 51049, 51128, 51179, 51234, 51336, 51345, 51641, 51649, 60589, 60627, 65007, 65009, 70727, 71702, 71705, 98772]}, 'taxonomy': {'entries': [{'_id': 39947, 'name': 'Oryza sativa Japonica Group'}], 'ancestors': [1, 2759, 3193, 3398, 4447, 4479, 4527, 4530, 4734, 33090, 35493, 38820, 58023, 58024, 78536, 131221, 131567, 147367, 147380, 359160, 1437183, 1437197, 1648021]}, 'familyRoot': {'entries': [{'_id': 2759, 'name': 'Eukaryota'}], 'ancestors': [1, 131567]}}, 'homology': {'gene_tree': {'id': 'EPlGT00920000145990', 'root_taxon_id': 2759, 'root_taxon_name': 'Eukaryota', 'duplications': [2759, 3193], 'representative': {'model': {'id': 'Os01g0100100', 'taxon_id': 39947, 'description': 'RabGAP/TBC domain containing protein. (Os01t0100100-01)'}}}, 'homologous_genes': {'ortholog_one2one': ['ORUFI01G00040', 'BGIOSGA002569', 'OBART01G00010', 'KN541486.1_FG001', 'OPUNC01G00010', 'OB01G10010', 'LPERR01G00020', 'F775_06225', 'TraesCS3A01G087000', 'BRADI_1g74660v3', 'TraesCS3B01G102400', 'TraesCS3D01G087200', 'HORVU3Hr1G015620', 'SETIT_000493mg', 'Zm00001d040234', 'SORBI_3003G110000', 'DCAR_029914', 'Dr08979', 'Solyc03g058470.1', 'PHAVU_005G155000g', 'PRUPE_1G200000', 'TCM_019958', 'Tp57577_TGAC_v2_gene6502', 'LR48_Vigan05g203300', 'CCACVL1_07405', 'Csa_4G010920', 'HannXRQ_Chr03g0089131', 'MTR_2g099185', 'AMTR_s00066p00146980', 'BVRB_3g057880', 'PHYPA_000100'], 'syntenic_ortholog_one2one': ['ONIVA01G00100', 'OGLUM01G00020'], 'ortholog_one2many': ['Bra007460', 'POPTR_017G023200v3', 'GLYMA_13G346200', 'VIT_00s0226g00160', 'GSBRNA2T00002271001', 'MANES_11G126700', 'GSBRNA2T00132657001', 'B456_005G126100', 'POPTR_007G132900v3', 'Bo4g019830', 'GLYMA_15G028000', 'GSMUA_Achr9G08200_001', 'GSBRNA2T00021969001', 'GSBRNA2T00121831001', 'AT3G59570', 'Bo4g194370', 'B456_011G138800', 'GSMUA_Achr6G28710_001', 'VIT_00s0226g00020', 'GSBRNA2T00155637001', 'scaffold_403239.1', 'AT2G43490', 'Bra037712', 'Bo8g094330', 'MANES_04G048100', 'Bra004766', 'fgenesh2_kg.5__2415__AT3G59570.1', 'GSBRNA2T00045107001', 'SELMODRAFT_416728', 'SELMODRAFT_442016', 'OSTLU_26794'], 'within_species_paralog': ['Os11g0587500', 'Os07g0525400', 'Os09g0528800', 'Os10g0518100', 'Os08g0547200', 'Os02g0709800', 'Os01g0908100', 'Os09g0515800']}}, 'bins': {'fixed_100': 3106, 'fixed_200': 6206, 'fixed_500': 15506, 'fixed_1000': 31008, 'uniform_1Mb': 31350, 'uniform_2Mb': 15789, 'uniform_5Mb': 6466, 'uniform_10Mb': 3342}, 'species_idx': 5},
...
'rapdb': {'ID': 'Os01g0100100', 'Description': 'RabGAP/TBC domain containing protein. (Os01t0100100-01)', 'Position': 'chr01:2983..10815', 'RAP-DB Gene Symbol Synonym(s)': '', 'RAP-DB Gene Name Synonym(s)': '', 'CGSNL Gene Symbol': '', 'CGSNL Gene Name': '', 'Oryzabase Gene Symbol Synonym(s)': '', 'Oryzabase Gene Name Synonym(s)': ''}, 'ic4r': {'All': '', 'OS Gene ID': 'Os01g0100100', 'LOC Gene ID': 'LOC_Os01g01010', 'Symbol': '-', 'Location': 'Chr1:2903-10817 (+)', 'Description': 'TBC domain containing protein, expressed', 'Box Plot': ''}}}

Process finished with exit code 0

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
* SNP-Seek
* funricegene (3 databse)
* MSU

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
