# -*- coding: utf-8 -*-
from MultiQuery import MultiQuery
import time
import re
import csv
import pickle

if __name__ == "__main__":
    test = MultiQuery()
    chromosome = ["chr01", "chr02", "chr03", "chr04", "chr05", "chr06", "chr07", "chr08", "chr09", "chr10", "chr11",
                  "chr12"]
    start = ["1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"]
    end = ["43270923", "35937250", "36413819", "35502694", "29958434", "31248787", "29697621", "28443022", "23012720"
        , "23207287", "29021106", "27531856"]
    iric_dict = dict()
    for i in range(len(chromosome)):
        t = time.time()
        file_id = test.search_gene(chro=chromosome[i], start_pos=start[i],
                                   end_pos=end[i], dbs=["iric"], save_path=None)
        iric_dict.update(file_id)
        print("Time for search gene:", time.time() - t)
    id_dict = dict()
    loc_dict = dict()
    for iric_name in iric_dict.keys():
        if iric_dict[iric_name]["msu7Name"] != None:
            for loc in iric_dict[iric_name]["msu7Name"]:
                if loc not in loc_dict.keys():
                    loc_dict.setdefault(loc, iric_name)
        if iric_dict[iric_name]["raprepName"] != None:
            for ids in iric_dict[iric_name]["raprepName"]:
                if ids not in id_dict.keys():
                    id_dict.setdefault(ids, iric_name)

    with open('./support/iric_dict', 'wb') as f:
        pickle.dump(iric_dict, f)
        f.close()
    with open('./support/loc_dict', 'wb') as f:
        pickle.dump(loc_dict, f)
        f.close()
    with open('./support/id_dict', 'wb') as f:
        pickle.dump(id_dict, f)
        f.close()
