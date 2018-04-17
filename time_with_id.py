import pandas as pd
import re
from urllib import urlopen
from bs4 import BeautifulSoup
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False
html = urlopen("http://www.lipidmaps.org/data/standards/fa_stds5_LC.php").read()
soup = BeautifulSoup(html)
table = soup.find("table", attrs={"class":"datatable"})
a = str(table)
y = 0
di = {}
ret = []
for x in re.split(u"(\|)|\-|\:|\;|\?|\(|\)|\*|\=|\\|\&|\/|\<|\>|\[|\]|\{|\}|\#|\+|\&|\%20|\_|\&nbsp|(\')|(\")", a):
    if x != None:
        if x[:4] == "LMFA":
            di[x] = 1
        elif is_number(x):
            if float(x) % 1 != 0:
                y += 1
                ret.append(x) 
print len(di)
a = 0
for key in di:
    print key, ret[a], ret[a+1]
    a += 2