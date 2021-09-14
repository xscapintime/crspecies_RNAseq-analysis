# %%
import os,glob
import sys
import pandas as pd
from bs4 import BeautifulSoup

# %%
path = '*_MotifOutput'
file = glob.glob(os.path.join(path,'homerResults.html'))

# %%
# for getting the header from
# the HTML file
# https://www.geeksforgeeks.org/convert-html-table-into-csv-file-in-python/
for f in file:
    data = []
    list_header = []
    soup = BeautifulSoup(open(f),'html.parser')
    header = soup.find_all("table")[0].find("tr")

    for items in header:
        try:
            list_header.append(items.get_text())
        except:
            continue

    # for getting the data 
    HTML_data = soup.find_all("table")[0].find_all("tr")[1:]
    
    for element in HTML_data:
        sub_data = []
        for sub_element in element:
            try:
                sub_data.append(sub_element.get_text())
            except:
                continue
        data.append(sub_data)

    dataFrame = pd.DataFrame(data = data, columns = list_header)
    dataFrame.replace({'\n': '', 'More Information | Similar Motifs Found':' ', r'\|':''}, regex=True,inplace=True)
    dataFrame[['Best Match/Details', 'score']] = dataFrame['Best Match/Details'].str.replace(' ','').str.split(r'\(0.', expand = True)
    dataFrame['score'] = '0.' + dataFrame['score'].str.replace(r'\)', '')

    # export
    file_name = os.path.split(f)[1].split('.html')[0] + '.txt'
    dataFrame.to_csv(os.path.split(f)[0] + '/' + file_name, sep='\t', header=True, index=False)