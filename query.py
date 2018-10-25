import helper
from bs4 import BeautifulSoup
import json
import regex
import pandas as pd
import csv
import requests

def fetch_description(db):
    # Fetch database description
    database_description = BeautifulSoup(open(
        "database-description.xml").read(), "xml").findAll("database", dbname=db.lower())
    if not database_description:
        raise ValueError('Database Not Found')
    return database_description

# Compose and execute query to database
def execute_query(db, qfields=[], verbose=False):
    # Get query qfields list
    fields = db[0].findAll("field")
    # Prepare URL
    link = db[0].findAll("link")[0]["stern"]
    # Compile URL
    if link[:4]=='http':
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
                    db[0].findAll("link")[0]["aft"]
                if verbose: print(link)
            return helper.connectionError(link)
    else:
        return open(link)

def query(db, qfields=[], outputFormat="dict", outputFile=None, verbose=False):
    
    # Fetch database description
    database_description = fetch_description(db)

    # Get Headers list
    headers = []
    for header in database_description[0].findAll("header"):
        headers.append(header.text)

    res = execute_query(database_description, qfields, verbose)
    if verbose:
        print(res.content)
    
    # Handle HTML based query
    if(database_description[0]["type"] == "text/html"):
        # Handling Connection
        ret = BeautifulSoup(res.content, "lxml")
        if verbose: 
            print("return value:", ret)
            print(database_description[0].findAll("data_struct")[0]["identifier"])
        if database_description[0].findAll("data_struct")[0]["identifier"] != "": 
            data = ret.findAll(database_description[0].findAll("data_struct")[0]["indicator"],
                        {database_description[0].findAll("data_struct")[0]["identifier"]:
                        database_description[0].findAll("data_struct")[0]["identification_string"]})
        else: data = ret.findAll(database_description[0].findAll("data_struct")[0]["indicator"])
        result = []
        if data != []:
            reg = regex.compile(database_description[0].findAll(
                "prettify")[0].text, regex.MULTILINE)
            replaceBy = database_description[0].findAll(
                "prettify")[0]['replaceBy']
            for dataLine in data[0].findAll(database_description[0].findAll("data_struct")[0]["line_separator"]):
                dict = {}
                i = 0
                for dataCell in dataLine.findAll(database_description[0].findAll("data_struct")[0]["cell_separator"]):
                    dataFormat = regex.sub(replaceBy+'+',' ',dataCell.text)
                    if(i<headers.__len__()):
                        dict[headers[i]] = dataFormat
                    i += 1
                if dict == {}:
                    continue
                dict.pop("", None)
                if verbose: print(dict)
                result.append(dict)
    # Handle JSON based query
    elif(database_description[0]["type"] == "text/JSON"):
        # Return as a List of Dictionary
        result = json.loads(res.content.decode("UTF-8"))
        if verbose: print (result)
    
    # Handle csv based DB
    # Auto detect and support local csv databases
    fields = database_description[0].findAll("field")
    if(database_description[0]["type"] == "text/csv"):
        if type(res) is requests.models.Response:
            ret = csv.reader(res.content.decode(database_description[0]["encoding"]).splitlines(), delimiter=list(database_description[0]["deli"])[0], quoting=csv.QUOTE_NONE)
        else:
            ret = csv.reader(res, delimiter=list(database_description[0]["deli"])[0], quoting=csv.QUOTE_NONE)
        result = []
        for row in ret:
            i = 0
            dict = {}
            for header in headers:
                dict[header] = row[i]
                i += 1
            f = 0
            for field in fields:
                match = True
                if (dict[field.text] != qfields[f]) & (qfields[f] != ""):
                    match = False
                    break
                f += 1
            if match: result.append(dict)
    
    # Handle different Output format
    df = pd.DataFrame(pd.io.json.json_normalize(result))
    if(outputFormat == "dict"):
        return result
    elif(outputFormat == "pandas"):
        return df
    elif(outputFormat == "json"):
        if (outputFile != None):
            print("Query exported to ", outputFile)
        return df.to_json(outputFile)
    elif(outputFile == None):
        print("Please specify a destination")
        return []
    if(outputFormat == "csv"):
        df.to_csv(outputFile)
        print("Query exported to ", outputFile)
        return []
    if(outputFormat == "excel"):
        df.to_excel(outputFile)
        print("Query exported to ", outputFile)
        return []