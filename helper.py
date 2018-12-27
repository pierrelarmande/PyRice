"""Helper module
Handle supportive function such as connection handling, path checking and error handling
"""
# -*- coding: utf-8 -*-
import os
import requests
import sys
import time
from bs4 import BeautifulSoup

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', elasped_time = None):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    if elasped_time != None:
        now = time.time() - elasped_time
        now_str = "Elasped time: %02d:%02d:%02d" % (now // 3600, (now % 3600 // 60), (now % 60 // 1))
        print('\r%s |%s| %s%% %s, %s' % (prefix, bar, percent, suffix, now_str), end = '\r')
    else:
        now_str = ""
        print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()

def existFile(pathToFile):
    """
    :param pathToFile: entire path to the file
    :return: return True if the file already exist, else return False
    """
    return (os.path.isfile(pathToFile))


def formatPathToFile(nameFile):
    """
    :param nameFile: name of the file with its extension
    :return: return the entire path to the file
    """

    # remove char until '/'
    pathToFile = os.path.dirname(__file__)
    #while not (pathToFile.endswith('/')):
    #    pathToFile = pathToFile[0:-1]

    pathToFile += '/resources/'+nameFile
    return pathToFile


def loadFileURL(url, localFilename=''):
    """
    Download the file located in the rapdb download page

    :param nameFile: name of the file (the all path to the file if you want to save the file in another folder)
    :param url: url where is located the file

    """

    # Fetch the file by the url and decompress it
    r = requests.get(url, stream=True)
    total = r.headers['content-length']
    print(total)
    progress = 0
    localFilename = url.split('/')[-1] and localFilename == ""
    with open(localFilename, 'wb') as f:
        for chunk in r.iter_content(chunk_size = 1024):
            progress+=len(chunk)
            # print(progress/int(total)*100)
            printProgressBar(progress, total, prefix='Downloading', suffix='', decimals=2)
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)
        f.close()
    return 1

def connectionError(link, data=""):

    """
    Return requests.get(link) with post request and test website issues

    :param link: URL
    :param data: data to give to the form
    :return: requests.get(link)
    """
    try:
        headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.0; WOW64; rv:24.0) Gecko/20100101 Firefox/24.0'}
        if data!= "":
            res = requests.post(link, data=data, headers=headers)
        else:
            res = requests.get(link, allow_redirects=False,timeout=5)
        if res.status_code != 200:
            print('Server Error: ' + str(res.status_code) + '\n' + 'For url:' + link)
            #raise Exception('Server Error: ' + str(res.status_code) + '\n' + 'For url:' + link)
            # sys.exit(1)
        return res
    except requests.exceptions.RequestException:
        print("Can't connect:",link)
        return
        #raise Exception("Internet Connection error")
        # sys.exit(1)