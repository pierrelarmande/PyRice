from multiprocessing import Process,Lock,Manager,active_children
from bs4 import BeautifulSoup
import helper
import json
import regex
import urllib
import requests
import time
import pickle
import csv
import pandas as pd
from pathlib import Path
import os
import copy


class MultiQuery():
    def __init__(self):
        self.result = None
        with open('./support/loc_dict', 'rb') as f:
            self.loc_dict = pickle.load(f)
            f.close()
        with open('./support/id_dict', 'rb') as f:
            self.id_dict = pickle.load(f)
            f.close()

    def fetch_description(self, db):
        # Fetch database description
        database_description = BeautifulSoup(open(
            "database-description.xml").read(), "xml").findAll("database", dbname=db.lower())
        if not database_description:
            raise ValueError('Database Not Found')
        return database_description

    def query(self, result, db, qfields=[], outputFormat="dict", outputFile=None, verbose=False):
        #         self.lock.acquire()
        # Fetch database description
        database_description = self.fetch_description(db)
        # Get Headers list
        headers = []
        for header in database_description[0].find_all("header"):
            headers.append(header.text)

        res = self.execute_query(database_description, qfields, verbose)
        if res.status_code != 200:
            return
        # Handle HTML based query
        if (database_description[0]["type"] == "text/html"):
            # Handling Connection
            ret = BeautifulSoup(res.content, "lxml")
            if verbose:
                print("return value:", ret)
                print(database_description[0].findAll("data_struct")[0]["identifier"])
            if database_description[0].findAll("data_struct")[0]["identifier"] != "":
                data = ret.find_all(database_description[0].find_all("data_struct")[0]["indicator"],
                                    {database_description[0].find_all("data_struct")[0]["identifier"]:
                                    database_description[0].find_all("data_struct")[0]["identification_string"]})
            else:
                data = ret.find_all(database_description[0].findAll("data_struct")[0]["indicator"])
                # result = []
            if data != []:
                reg = regex.compile(database_description[0].findAll(
                    "prettify")[0].text, regex.MULTILINE)
                replaceBy = database_description[0].findAll(
                    "prettify")[0]['replaceBy']
                for dataLine in data[0].findAll(database_description[0].findAll("data_struct")[0]["line_separator"]):
                    dict_ = {}
                    i = 0
                    for dataCell in dataLine.findAll(
                            database_description[0].findAll("data_struct")[0]["cell_separator"]):
                        dataFormat = regex.sub(replaceBy + '+', ' ', dataCell.text)
                        if (i < len(headers)):
                            dict_[headers[i]] = dataFormat
                        i += 1
                    if dict_ == {}:
                        continue
                    dict_.pop("", None)
                    if verbose:
                        print(dict_)
                    self.result[db].setdefault(qfields[-1], dict_)
                    #self.result[db].append([qfields[-1], dict_])
        # Handle JSON based query
        elif (database_description[0]["type"] == "text/JSON"):
            # Return as a List of Dictionary
            if len(res.content) > 10:
                self.result[db].setdefault(qfields[-1],json.loads(res.content))
            #self.result[db].append([qfields[-1], json.loads(res.content)])
            if verbose: print(self.result[db])
        # Handle csv based DB
        # Auto detect and support local csv databases
        fields = database_description[0].findAll("field")
        if (database_description[0]["type"] == "text/csv"):
            if type(res) is requests.models.Response:
                ret = csv.reader(res.content.decode(database_description[0]["encoding"]).splitlines(),
                                 delimiter=list(database_description[0]["deli"])[0], quoting=csv.QUOTE_NONE)
            else:
                ret = csv.reader(res, delimiter=list(database_description[0]["deli"])[0], quoting=csv.QUOTE_NONE)
            l_res = set()
            tmp = []
            for row in ret:
                tup_row = tuple(row)
                if tup_row not in l_res:
                    l_res.add(tup_row)
                    i = 0
                    dict_ = {}
                    for header in headers:
                        dict_[header] = row[i]
                        i += 1
                    f = 0
                    for field in fields:
                        match = True
                        if (dict_[field.text] != qfields[f]) & (qfields[f] != ""):
                            match = False
                            break
                        f += 1
                    if match:
                        tmp.append(dict_)
            if len(tmp)>0:
                self.result[db].setdefault((qfields[0],qfields[1]),tmp)
            # for t in tmp:
                #self.result[db].append([(qfields[0], qfields[1]), t])
        return self.result

    def save_file(self,folder_path):
        data_folder = folder_path+"data/"
        gene_folder = folder_path+"gene/"
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)
        if not os.path.exists(gene_folder):
            os.makedirs(gene_folder)
        #ouput database
        test = copy.deepcopy(self.result)
        filter = ["db_type", "gene_idx", "location", "xrefs", "gene_structure", "bins", "annotations", "homology"]
        for i in test.keys():
            for k, v in test[i].items():
                if k == "Gramene":
                    for att in filter:
                        v[0].pop(att, None)
        my_db = data_folder +"db" + '.csv'
        df = pd.DataFrame.from_dict(test, orient='index')
        with open(my_db, 'w') as f:
            df.to_csv(f, header=True)
            f.close()
        #output gene
        test = copy.deepcopy(self.result)
        filter = ["gene_structure","bins",]
        annotations = ["taxonomy","familyRoot"]
        homology = ["gene_tree"]
        homologous_genes =["syntenic_ortholog_one2one","ortholog_one2many"]
        for i in test.keys():
            for k,v in test[i].items():
                   if k == "Gramene":
                        for att in filter:
                            v[0].pop(att,None)
                        if "annotations" in v[0].keys():
                            for att in annotations:
                                v[0]["annotations"].pop(att,None)
                        if "homology" in v[0].keys():
                            for att in homology:
                                v[0]["homology"].pop(att,None)
                            if "homologous_genes" in v[0].keys():
                                for att in homologous_genes:
                                    v[0]["homology"]["homologous_genes"].pop(att,None)
        for key,value in test.items():
            my_gene = gene_folder +key +'.txt'
            with open(my_gene,'w') as f:
                f.write(json.dumps(test[key],indent="\t"))
                f.close()

    def execute_query(self, db, qfields=[], verbose=False):
        #Get query qfields list
        fields = db[0].find_all("field")
        # Prepare URL
        link = db[0].find_all("link")[0]["stern"]
        # Compile URL
        if link[:4] == 'http':
            if db[0]["method"] == "POST":
                i = 0
                for field in fields:
                    data = {field.text: qfields[i]}
                    i += 1
                return helper.connectionError(link, data)
            elif db[0]["method"] == "GET":
                query_string = ""
                if db[0]["type"] != "text/csv":
                    i = 0
                    for field in fields:
                        # Detect controller field (always first field)
                        if "lowercase" in field:
                            print(qfields[i].lower())
                        if field.text == "":
                            query_string += qfields[i] + "?"
                        # All other fields are query fields
                        else:
                            query_string += field.text + field["op"] + qfields[i] + "&"
                        i += 1
                    query_string = query_string[:-1]
                    link += query_string + \
                            db[0].find_all("link")[0]["aft"]
                    if verbose: print(link)
                return helper.connectionError(link)
        else:
            return open(link)

    def query_all(self, idents, locs, dbs='all'):
        set_ids,set_locs,idents, locs = self.check_gene(idents, locs)
        manager = Manager()
        self.result = manager.dict()
        if dbs == 'all':
            name_db = ["oryzabase", "Gramene", "funricegene_genekeywords",
                       "funricegene_faminfo", "msu", "rapdb","ic4r",
                       "funricegene_geneinfo"]
        else:
            name_db = dbs
        number_process = 100
        p =[None for i in range(number_process)]
        count=0
        for db in name_db:
            if db not in self.result.keys():
                #self.result.setdefault(db, manager.list())
                self.result.setdefault(db, manager.dict())
            if db == "rapdb" or db == 'oryzabase' or db == "Gramene" or db =="ic4r":
                for ident in set_ids:
                    p[count] = Process(target=self.query, args=(self.result, db, [ident],))
                    count+=1
                    p[count-1].start()
                #     list_proccess.append(p)
                # for p in list_proccess:
                #     p.join();
            elif db == "msu":
                for loc in set_locs:
                    p[count] = Process(target=self.query, args=(self.result, db, [loc],))
                    count+=1
                    p[count-1].start()
                #     list_proccess.append(p)
                # for p in list_proccess:
                #     p.join();
            elif db == "funricegene_genekeywords" or db == "funricegene_faminfo" or db =="funricegene_geneinfo":
                for ident in idents.keys():
                    for loc in idents[ident]:
                        p[count] = Process(target=self.query, args=(self.result, db, [ident, loc],))
                        count+=1
                        p[count-1].start()
                    #     list_proccess.append(p)
                    # for p in list_proccess:
                    #     p.join();
                for loc in locs.keys():
                    for ident in locs[loc]:
                        p[count] = Process(target=self.query, args=(self.result, db, [ident, loc],))
                        count+=1
                        p[count-1].start()
                    #     list_proccess.append(p)
        for i in range(0,count):
            p[i].join()
        new_result = dict()
        new_name_db=[]
        fun_db =[]
        if dbs == 'all':
            new_name_db = ["oryzabase", "Gramene", "rapdb", "ic4r"]
            fun_db = ["funricegene_genekeywords", "funricegene_faminfo", "funricegene_geneinfo"]
        else:
            for db in name_db:
                if db == "rapdb" or db == 'oryzabase' or db == "Gramene" or db == "ic4r":
                    new_name_db.append(db)
                else:
                    fun_db.append(db)
        for ident in idents.keys():
            id_data = dict()
            for db in new_name_db:
                #check ident in {db:{id:data}}
                if ident in self.result[db].keys():
                    #get data
                    id_data.setdefault(db,self.result[db][ident])
            #create dict with new_key : id-loc
            new_key = ident
            loc_data = dict()
            #get locs for create new_key
            for loc in idents[ident]:
                # search id-loc in fun
                new_key += "-" + loc
                # search loc in msu
                if "msu" in name_db:
                    if loc in self.result["msu"].keys():
                        loc_data.setdefault("msu",self.result["msu"][loc])
                        id_data.update(loc_data)
                fun_data = dict()
                for fun in fun_db:
                    if (ident,loc) in self.result[fun].keys():
                        fun_data.setdefault(fun,self.result[fun][ident,loc])
                        id_data.update(fun_data)
            new_result.setdefault(new_key,id_data)
        for loc in locs.keys():
            if loc in self.result["msu"].keys():
                new_result.setdefault(loc,self.result["msu"][loc])
        #print("new_result",new_result)
        with open('./support/data.json', 'w') as fp:
            json.dump(new_result, fp)
        self.result = new_result
        #self.save_file("./data/")
        return self.result

    def check_gene(self, idents, locs):
        true_ids = dict()
        true_locs = dict()
        set_locs = set() #check locs
        set_ids = set() #check_ids
        for ident in idents:
            if ident in self.id_dict.keys():
                set_ids.add(ident)
                true_ids.setdefault(ident, self.id_dict[ident])
                for loc in self.id_dict[ident]:
                    set_locs.add(loc)
            #print("{} : {}".format(ident,self.id_dict[ident]))
            else:
                print("Can't find Id: ", ident)
        for loc in locs:
            if loc not in set_locs:
                if loc in self.loc_dict.keys():
                    true_locs.setdefault(loc, set())
                    for ident in self.loc_dict[loc]:
                        if ident not in set_ids:
                            if ident == 'None':
                                true_locs[loc].add(ident)
                            else:
                                set_ids.add(ident)
                                true_ids.setdefault(ident,self.id_dict[ident])
                                true_locs.pop(loc,None)
                    #print("{}:{}".format(loc,self.loc_dict[loc]))
                else:
                    print("Can't find Loc: ", loc)
        #print(set_locs,set_ids)
        #print(true_ids,true_locs)
        return set_ids,set_locs,true_ids, true_locs

    def search_gene(self, chro, start_pos, end_pos, dbs='all'):
        manager = Manager()
        self.result = manager.dict()
        if dbs == 'all':
            name_db = ["rap", "msu7", "iric"]
        else:
            name_db = dbs
        #self.result.setdefault('snpseek', manager.list())
        self.result.setdefault('snpseek', manager.dict())
        number_process = len(name_db)
        list_process = [None for i in range(number_process)]
        for i in range(number_process):
            list_process[i] = Process(target=self.query,
                        args=(self.result, "snpseek", [str(chro), str(start_pos), str(end_pos), name_db[i] ],))
            list_process[i].start()
        for process in list_process:
            process.join()
        return self.result

    def f(self, result, db, qfields=[], outputFormat="dict", outputFile=None, verbose=False):
        self.lock.acquire()
        self.result.append({'db': db, 'gene': qfields})
        self.lock.release()
        return self.result

#if __name__ == "__main__":
    # ids, locs = MultiQuery().check_gene(idents=["Os08g0164400", "Os07g0586200"]
    #                                     ,locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914"])
    # print(ids, locs)

    # t = time.time()
    # multi_query = MultiQuery()
    # test = multi_query.query_all(idents=["Os08g0164400", "Os07g0586200","Os01g0100900"]
    #                                          ,locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914","LOC_Os01g01019"])
    #test = multi_query.query_all(idents=[],locs=["LOC_Os01g01019"],dbs=["Gramene","msu"])
    # multi_query.save_file("./result/")
    # print(time.time() - t)
    # filter = ["db_type","gene_idx","location","xrefs","gene_structure","bins","annotations","homology"]
    # for i in test.keys():
    #     print("Gene:", i)
    #     for k,v in test[i].items():
    #         if k == "Gramene":
    #             for att in filter:
    #                 v[0].pop(att,None)
    #         print(k,v)
    #     print("\n------------------------")
    # filter = ["gene_structure","bins",]
    # annotations = ["taxonomy","familyRoot"]
    # homology = ["gene_tree"]
    # homologous_genes =["syntenic_ortholog_one2one","ortholog_one2many"]
    # for i in test.keys():
    #     print("Gene:", i)
    #     for k,v in test[i].items():
    #            if k == "Gramene":
    #                 for att in filter:
    #                     v[0].pop(att,None)
    #                 for att in annotations:
    #                     v[0]["annotations"].pop(att,None)
    #                 for att in homology:
    #                     v[0]["homology"].pop(att,None)
    #                 for att in homologous_genes:
    #                     v[0]["homology"]["homologous_genes"].pop(att,None)
    #         print(k,v)
    #     print("\n------------------------")
    # t = time.time()
    # test = MultiQuery().search_gene(chro="chr01", start_pos="1",
    #                                 end_pos="20000", dbs=["msu7"])
    # print("Time for search gene:",time.time() - t)
    # for i in test.keys():
    #     print("Database:", i)
    #     for k,v in test[i].items():
    #         print(k,v)
    #     print("\n------------------------")
    # id1 = []
    # loc1 = []
    # id2 = []
    # loc2 = []
    # for key,value in test.items():
    #     for db,list_gene in value.items():
    #         if db == 'msu7':
    #             id_att = "raprepName"
    #             loc_att = "msu7Name"
    #         elif db == 'rap':
    #             id_att = "uniquename"
    #             loc_att = "msu7Name"
    #         else:
    #             id_att = "raprepName"
    #             loc_att = "msu7Name"
    #         for gene in list_gene:
    #             if gene[id_att] !=None:
    #                 id1.append(gene[id_att])
    #                 #print(gene[id_att])
    #             if gene[loc_att] != None:
    #                 loc1.append(gene[loc_att])
    #                 #print(gene[loc_att])
    # print(id1, loc1)
    # t = time.time()
    # # id1.append('Os09g0360900')
    # multi_query = MultiQuery()
    # test = multi_query.query_all(idents=id1 ,locs=loc1)
    # multi_query.save_file("./result/")
    # print("Time for query:",time.time() - t)

