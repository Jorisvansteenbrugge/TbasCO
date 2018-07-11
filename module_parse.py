from HTMLParser import HTMLParser
import requests
import unicodedata

class myhtmlparser(HTMLParser):
    def __init__(self):
        self.reset()
        self.NEWTAGS = []
        self.NEWATTRS = []
        self.HTMLDATA = []
    def handle_starttag(self, tag, attrs):
        self.NEWTAGS.append(tag)
        self.NEWATTRS.append(attrs)
    def handle_data(self, data):
        self.HTMLDATA.append(data)
    def clean(self):
        self.NEWTAGS = []
        self.NEWATTRS = []
        self.HTMLDATA = []



def parse_HTML_list(html_list):
    out = []
    start = False
    for val in html_list:
        if 'Type' in val:
            start = False
        if start == True:
            out.append(val)
        elif 'Definition' in val:
            start = True

    # Pop empy lines
    clean_out = []
    for val in out:
        if '\n' in val:
            pass
        else:
            clean_out.append(val.encode('ascii'))

    return clean_out


def retrieve_module_versions(parsed_KOs):


def parseModule(module):
    module_specific = "http://www.genome.jp/kegg-bin/show_module?{}"
    r = requests.get(module_specific.format(module))
    r = r.text

    parser = myhtmlparser()
    parser.feed(r)
    # Extract data from parser

    data  = parser.HTMLDATA

    # Clean the parser
    parser.clean()

    parsed_KOs = parse_HTML_list(data)



parseModule('M00050')
