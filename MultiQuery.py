# -*- coding: utf-8 -*-
from multiprocessing import Process,Lock,Manager,active_children,Pool
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

    def query(self,iricname, db, qfields=[], outputFormat="dict", outputFile=None, verbose=False):
        """
        Query one gene with id or loc in each database
        :param iricname: iricname or id of gene
        :param db: name database in 9 databases
        :param qfields: list of loc,id
        :param outputFormat: None
        :param outputFile: None
        :param verbose: if True print for debug
        :return: None , support for query_iric and query_ids_locs
        """
        #         self.lock.acquire()
        # Fetch database description
        database_description = helper.fetch_description(db)
        # Get Headers list
        headers = []
        for header in database_description[0].find_all("header"):
            headers.append(header.text)

        res = helper.execute_query(database_description, qfields, verbose)
        if res == None:
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
                    #self.result[db].setdefault(qfields[-1], dict_)
                    if len(dict_) >0:
                        self.result[iricname].setdefault(db,dict_)
                    #self.result[db].append([qfields[-1], dict_])
        # Handle JSON based query
        elif (database_description[0]["type"] == "text/JSON"):
            # Return as a List of Dictionary
            if len(res.content) > 10:
                # self.result[db].setdefault(qfields[-1],json.loads(res.content.decode('utf-8')))
                if iricname == "snpseek":
                    self.result[iricname].setdefault(qfields[-1],json.loads(res.content.decode('utf-8')))
                else:
                    self.result[iricname].setdefault(db,json.loads(res.content.decode('utf-8'))[0])
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
            dict_2 = dict()
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
                        for k,v in dict_.items():
                            if k not in dict_2.keys():
                                dict_2.setdefault(k,set())
                                dict_2[k].add(v)
                            else:
                                dict_2[k].add(v)
                        tmp.append(dict_)
            if len(tmp)>0:
                #self.result[db].setdefault((qfields[0],qfields[1]),tmp)
                for k, v in dict_2.items():
                    dict_2[k] = repr(v);
                self.result[iricname].setdefault(db,dict_2)
            # for t in tmp:
                #self.result[db].append([(qfields[0], qfields[1]), t])

    #Do not use
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

    def query_iric(self,file_id,dbs='all',save_path =None):
        """
        Query with iricname
        :param file_id: result of function search_gene
        :param dbs: list databases (support 9 available databases)
        :param save_path: (str) path to save result after call function
        :return: a dictionary, format : gene:{database: attributes}
        """
        set_ids = set()
        set_locs = set()
        # if iricnames == 'all':
        #     for key,value in file.items():
        #         for l in value["msu7Name"]:
        #             set_locs.add(l)
        #         for i in value["raprepName"]:
        #             set_ids.add(i)
        # else:
        #     for iricname in iricnames:
        #         if iricname in file.keys():
        #             for l in file[iricname]["msu7Name"]:
        #                 set_locs.add(l)
        #             for i in file[iricname]["raprepName"]:
        #                 set_ids.add(i)
        manager = Manager()
        self.result = manager.dict()
        #Check support of database
        support_db = ["oryzabase", "Gramene", "funricegene_genekeywords",
                   "funricegene_faminfo", "msu", "rapdb", "ic4r",
                   "funricegene_geneinfo"]
        if dbs == 'all':
            name_db = support_db
        else:
            name_db=[]
            for db in dbs:
                if db not in support_db:
                    print("Don't support databse: ", db)
                else:
                    name_db.append(db)
        list_process = []
        #Query in multi_database
        i=1
        number_query = len(file_id)
        try:
            p = Pool(processes=8)
            for key, value in file_id.items():
                print("Query iricname: {} --- Step {}/{}".format(key, i, number_query))
                i+=1
                self.result.setdefault(key,manager.dict())
                for db in name_db:
                    if db == "rapdb" or db == 'oryzabase' or db == "Gramene" or db == "ic4r":
                        for ident in value["raprepName"]:
                            # p = Process(target=self.query, args=(key,db,[ident],))
                            # p.start()
                            # list_process.append(p)
                            p.apply_async(self.query,args=(key,db,[ident],))
                    elif db == "msu":
                        for loc in value["msu7Name"]:
                            # p = Process(target=self.query, args=(key,db, [loc],))
                            # p.start()
                            # list_process.append(p)
                            p.apply_async(self.query, args=(key,db, [loc],))
                    elif db == "funricegene_genekeywords" or db == "funricegene_faminfo" or db == "funricegene_geneinfo":
                        if len(value["raprepName"]) >0:
                            for ident in value["raprepName"]:
                                for loc in value["msu7Name"]:
                                    # p = Process(target=self.query, args=(key, db, [ident, loc],))
                                    # p.start()
                                    # list_process.append(p)
                                    p.apply_async(self.query, args=(key, db, [ident,loc],))
                        else:
                            for loc in value["msu7Name"]:
                                # p = Process(target=self.query, args=(key, db, ["", loc],))
                                # p.start()
                                # list_process.append(p)
                                p.apply_async(self.query, args=(key, db, ["",loc],))
            # for process in list_process:
            #     process.join()
        finally:
            p.close()
            p.join()
            #p.terminate()
        #change format to save file
        demo = dict()
        for key, value in self.result.items():
            demo.setdefault(key, dict())
            for subdb, data in value.items():
                demo[key].setdefault(subdb, data)
        if save_path != None:
            data_folder = save_path + "data/"
            gene_folder = save_path + "gene/"
            if not os.path.exists(data_folder):
                os.makedirs(data_folder)
            if not os.path.exists(gene_folder):
                os.makedirs(gene_folder)
            #ouput data/db.csv
            test = copy.deepcopy(demo)
            #filter allow attributes
            filter_ = ["db_type", "gene_idx", "location", "xrefs", "gene_structure", "bins", "annotations", "homology"]
            for i in test.keys():
                for k, v in test[i].items():
                    if k == "Gramene":
                        for att in filter_:
                            v.pop(att, None)
            my_db = data_folder + "db" + '.csv'
            #Convert output db.csv
            total_dict = dict()
            for iricname,databases in test.items():
                new_dict = dict()
                for db,data in databases.items():
                    new_data = dict()
                    for att,value in data.items():
                        if value != "":
                            new_data.setdefault(db+"."+att,value)
                    new_dict.update(new_data)
                total_dict.setdefault(iricname, new_dict)
            df = pd.DataFrame.from_dict(total_dict,orient='index')
            with open(my_db, 'w') as f:
                df.to_csv(f, header=True)
                f.close()
            # output gene/iricname.txt
            test = copy.deepcopy(demo)
            # filter allow attributes
            filter_ = ["gene_structure", "bins", ]
            annotations = ["taxonomy", "familyRoot"]
            homology = ["gene_tree"]
            homologous_genes = ["syntenic_ortholog_one2one", "ortholog_one2many"]
            for i in test.keys():
                for k, v in test[i].items():
                    if k == "Gramene":
                        for att in filter_:
                            v.pop(att, None)
                        if "annotations" in v.keys():
                            for att in annotations:
                                v["annotations"].pop(att, None)
                        if "homology" in v.keys():
                            for att in homology:
                                v["homology"].pop(att, None)
                            if "homologous_genes" in v.keys():
                                for att in homologous_genes:
                                    v["homology"]["homologous_genes"].pop(att, None)
            for iricname,databases in test.items():
                my_gene = gene_folder + iricname + '.txt'
                with open(my_gene, 'w') as f:
                    f.write(json.dumps(test[iricname], indent="\t"))
                    f.close()
        self.result = demo
        return self.result

    def query_ids_locs(self,idents, locs, dbs='all',save_path = None):
        """
        Query with id and loc of gene
        :param idents: list id of gene
        :param locs: list loc of gene
        :param dbs: list databases (support 9 available databases)
        :param save_path: (str) path to save result after call function
        :return: a dictionary, format : gene:{database: attribute}
        """
        set_ids,set_locs,idents,locs = self.check_gene(idents, locs)
        manager = Manager()
        self.result = manager.dict()
        support_db = ["oryzabase", "Gramene", "funricegene_genekeywords",
                   "funricegene_faminfo", "msu", "rapdb", "ic4r",
                   "funricegene_geneinfo"]
        if dbs == 'all':
            name_db = support_db
        else:
            name_db=[]
            for db in dbs:
                if db not in support_db:
                    print("Don't support databse: ", db)
                else:
                    name_db.append(db)
        list_process = []
        i = 1;
        number_query = len(idents)
        try:
            p = Pool(processes=8)
            for key, value in idents.items():
                print("Query ID: {} --- Step {}/{}".format(key,i,number_query))
                i+=1
                self.result.setdefault(key,manager.dict())
                for db in name_db:
                    if db == "rapdb" or db == 'oryzabase' or db == "Gramene" or db == "ic4r":
                        #p = Process(target=self.query, args=(self.result,key,db,[key],))
                        #p.start()
                        #list_process.append(p)
                        p.apply_async(self.query, args=(key,db,[key],))
                    elif db == "msu":
                        for loc in value:
                            #p = Process(target=self.query, args=(self.result,key,db, [loc],))
                            #p.start()
                            #list_process.append(p)
                            p.apply_async(self.query, args=(key,db, [loc],))
                    elif db == "funricegene_genekeywords" or db == "funricegene_faminfo" or db == "funricegene_geneinfo":
                        if len(value) >0:
                            for loc in value:
                                #p = Process(target=self.query, args=(self.result,key, db, [key,loc],))
                                #p.start()
                                #list_process.append(p)
                                p.apply_async(self.query, args=(key, db, [key,loc],))
                        else:
                            for loc in value:
                                #p = Process(target=self.query, args=(self.result,key, db, ["", loc],))
                                #p.start()
                                #list_process.append(p)
                                p.apply_async(self.query, args=(key, db, [key,loc],))
            if "msu" in name_db:
                print("Query LOC with out ID")
                i=1;
                number_query = len(locs)
                for loc in locs.keys():
                    self.result.setdefault(loc,manager.dict())
                    print("Query LOC: {} --- Step {}/{} ".format(loc,i,number_query))
                    i+=1;
                    #p = Process(target=self.query, args=(self.result,loc, "msu", [loc],))
                    #p.start()
                    #list_process.append(p)
                    p.apply_async(self.query, args=(loc, "msu", [loc],))
        finally:
            p.close()
            p.join()
        # for process in list_process:
        #     process.join()
        # change format to save file
        demo = dict()
        for key, value in self.result.items():
            demo.setdefault(key, dict())
            for subdb, data in value.items():
                demo[key].setdefault(subdb, data)
        if save_path != None:
            data_folder = save_path + "data/"
            gene_folder = save_path + "gene/"
            if not os.path.exists(data_folder):
                os.makedirs(data_folder)
            if not os.path.exists(gene_folder):
                os.makedirs(gene_folder)
            # ouput data/db.csv
            test = copy.deepcopy(demo)
            # filter allow attributes
            filter_ = ["db_type", "gene_idx", "location", "xrefs", "gene_structure", "bins", "annotations",
                       "homology"]
            for i in test.keys():
                for k, v in test[i].items():
                    if k == "Gramene":
                        for att in filter_:
                            v.pop(att, None)
            my_db = data_folder + "db" + '.csv'
            # Convert output db.csv
            total_dict = dict()
            for iricname, databases in test.items():
                new_dict = dict()
                for db, data in databases.items():
                    new_data = dict()
                    for att, value in data.items():
                        if value != "":
                            new_data.setdefault(db + "." + att, value)
                    new_dict.update(new_data)
                total_dict.setdefault(iricname, new_dict)
            df = pd.DataFrame.from_dict(total_dict, orient='index')
            with open(my_db, 'w') as f:
                df.to_csv(f, header=True)
                f.close()
            # output gene/iricname.txt
            test = copy.deepcopy(demo)
            # filter allow attributes
            filter_ = ["gene_structure", "bins", ]
            annotations = ["taxonomy", "familyRoot"]
            homology = ["gene_tree"]
            homologous_genes = ["syntenic_ortholog_one2one", "ortholog_one2many"]
            for i in test.keys():
                for k, v in test[i].items():
                    if k == "Gramene":
                        for att in filter_:
                            v.pop(att, None)
                        if "annotations" in v.keys():
                            for att in annotations:
                                v["annotations"].pop(att, None)
                        if "homology" in v.keys():
                            for att in homology:
                                v["homology"].pop(att, None)
                            if "homologous_genes" in v.keys():
                                for att in homologous_genes:
                                    v["homology"]["homologous_genes"].pop(att, None)
            for iricname, databases in test.items():
                my_gene = gene_folder + iricname + '.txt'
                with open(my_gene, 'w') as f:
                    f.write(json.dumps(test[iricname], indent="\t"))
                    f.close()
        self.result = demo
        return self.result
        # new_result = dict()
        # new_name_db=[]
        # fun_db =[]
        # if dbs == 'all':
        #     new_name_db = ["oryzabase", "Gramene", "rapdb", "ic4r"]
        #     fun_db = ["funricegene_genekeywords", "funricegene_faminfo", "funricegene_geneinfo"]
        # else:
        #     for db in name_db:
        #         if db == "rapdb" or db == 'oryzabase' or db == "Gramene" or db == "ic4r":
        #             new_name_db.append(db)
        #         else:
        #             fun_db.append(db)
        # for ident in idents.keys():
        #     id_data = dict()
        #     for db in new_name_db:
        #         #check ident in {db:{id:data}}
        #         if ident in self.result[db].keys():
        #             #get data
        #             id_data.setdefault(db,self.result[db][ident])
        #     #create dict with new_key : id-loc
        #     new_key = ident
        #     loc_data = dict()
        #     #get locs for create new_key
        #     for loc in idents[ident]:
        #         # search id-loc in fun
        #         new_key += "-" + loc
        #         # search loc in msu
        #         if "msu" in name_db:
        #             if loc in self.result["msu"].keys():
        #                 loc_data.setdefault("msu",self.result["msu"][loc])
        #                 id_data.update(loc_data)
        #         fun_data = dict()
        #         for fun in fun_db:
        #             if (ident,loc) in self.result[fun].keys():
        #                 fun_data.setdefault(fun,self.result[fun][ident,loc])
        #                 id_data.update(fun_data)
        #     new_result.setdefault(new_key,id_data)
        # for loc in locs.keys():
        #     if loc in self.result["msu"].keys():
        #         new_result.setdefault(loc,self.result["msu"][loc])
        # #print("new_result",new_result)
        # with open('./support/data.json', 'w') as fp:
        #     json.dump(new_result, fp)
        # self.result = new_result
        # #self.save_file("./data/")
        # return self.result

    def check_gene(self, idents, locs):
        """
        Check gene before query
        :param idents: list id of gene
        :param locs: list loc of gene
        :return:
            set_ids: set of id (existed)
            set_locs: set of locs (exitsted)
            true_ids: a dictionary, format: {id:loc}
            true_locs: a dictionary, format: {loc:None}
        """
        true_ids = dict()
        true_locs = dict()
        set_locs = set() #check locs
        set_ids = set() #check_ids
        if len(idents) > 1:
            for ident in idents:
                if ident in self.id_dict.keys():
                    set_ids.add(ident)
                    true_ids.setdefault(ident, self.id_dict[ident])
                    for loc in self.id_dict[ident]:
                        set_locs.add(loc)
                #print("{} : {}".format(ident,self.id_dict[ident]))
                else:
                    print("Can't find ID: ", ident)
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
                    print("Can't find LOC: ", loc)
        #print(set_locs,set_ids)
        #print(true_ids,true_locs)
        return set_ids,set_locs,true_ids, true_locs

    def search_gene(self, chro, start_pos, end_pos, dbs='all',save_path = None):
        """
        Search gene in snpseek
        :param chro: (str) chromosome
        :param start_pos: (str)
        :param end_pos: (str)
        :param dbs: list of database: ["rap", "msu7", "iric"]
        :param save_path: path to save result after call function
        :return: a dictionary, format: iricname:{{msu7Name:LOC_Os..},{raprepName:Os}}
        """
        manager = Manager()
        self.result = manager.dict()
        if dbs == 'all':
            name_db = ["rap", "msu7", "iric"]
        else:
            name_db = dbs
        self.result.setdefault('snpseek', manager.dict())
        list_process = []
        try:
            p = Pool(processes=8)
            for i in range(len(dbs)):
                # p = Process(target=self.query,
                #             args=("snpseek","snpseek",[str(chro), str(start_pos), str(end_pos), name_db[i]],))
                # p.start()
                # list_process.append(p)
                p.apply_async(self.query,args=("snpseek","snpseek",[str(chro), str(start_pos), str(end_pos), name_db[i]],))
        finally:
            p.close()
            p.join()
            # p.terminate()
        # for process in list_process:
        #     process.join()
        item_dict = dict()
        test = copy.deepcopy(self.result)
        # for i in test.keys():
        #     print("Gene:", i)
        #     for k,v in test[i].items():
        #         print(k,v)
        #     print("\n------------------------")
        key_att = "iricname"
        id_att = "raprepName"
        id_att_sup = "rappredName"
        loc_att = "msu7Name"
        #Build dictionary iricname:{{msu7Name:LOC_Os..},{raprepName:Os}}
        for key, value in test.items():
            for db, list_gene in value.items():
                for gene in list_gene:
                    if gene[key_att] not in item_dict.keys():
                        item_dict.setdefault(gene[key_att],dict())
                        item_dict[gene[key_att]].setdefault("msu7Name",set())
                        item_dict[gene[key_att]].setdefault("raprepName", set())
                    if gene[loc_att] != None:
                        #    print(gene[loc_att],type(gene[loc_att]))
                        for l in gene[loc_att].split(','):
                            item_dict[gene[key_att]]["msu7Name"].add(l)
                    if gene[id_att] != None:
                        #    print(gene[id_att],type(gene[id_att]))
                        for g in gene[id_att].split(','):
                            item_dict[gene[key_att]]["raprepName"].add(g)
                        continue
                    if gene[id_att_sup] != None:
                        #    print(gene[id_att],type(gene[id_att]))
                        for g in gene[id_att_sup].split(','):
                            item_dict[gene[key_att]]["raprepName"].add(g)
        #Save output
        if save_path != None:
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            my_db = save_path + "items.csv"
            df = pd.DataFrame.from_dict(item_dict, orient='index')
            with open(my_db, 'w') as f:
                df.to_csv(f, header=True)
                f.close()
        return item_dict

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

