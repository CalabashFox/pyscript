from bs4 import BeautifulSoup
import urllib
import urllib3
import json

search_url = 'http://cancer.sanger.ac.uk/cosmic/search/genes'

search_params = {
    'export': 'json',
    'sEcho': 1,
    'iColumns': 4,
    'sColumns': '',
    'iDisplayStart': 0,
    'iDisplayLength': 100,
    'sSearch': '',
    'bRegex': False,
    'iSortCol_0': 0,
    'sSortDir_0': 'asc',
    'iSortingCols': 1
}

fusion_url = 'http://cancer.sanger.ac.uk/cosmic/gene/fusions'

fusion_params = {
    'src': 'gene',
    'coords': 'AA:AA',
    'dr': '',
    'gd': '',
    'all_data': '',
    'start': 1,
    'export': 'json',
    'sEcho': 1,
    'iColumns': 6,
    'sColumns': '',
    'iDisplayStart': 0,
    'iDisplayLength': 600,
    'sSearch': '',
    'bRegex': False,
    'iSortCol_0': 0,
    'sSortDir_0': 'asc',
    'iSortingCols': 1
}

tissue_url = 'http://cancer.sanger.ac.uk/cosmic/gene/tissue'

tissue_params = {
    'src': 'gene',
    'coords': 'AA:AA',
    'dr': '',
    'gd': '',
    'all_data': '',
    'start': 1
}

frameshift_url = 'http://cancer.sanger.ac.uk/cosmic/gene/positive'

frameshift_params = {
    'src': 'gene',
    'coords': 'AA:AA',
    'gd': '',
    'all_data': '',
    'start': 1,
    'export': 'json',
    'sEcho': 1,
    'iColumns': 6,
    'sColumns': '',
    'iDisplayStart': 0,
    'iDisplayLength': 600,
    'sSearch': 'upper aerodigestive tract',
    'bRegex': False,
    'iSortCol_0': 0,
    'sSortDir_0': 'asc',
    'iSortingCols': 1
}

mut_overview_url = 'http://cancer.sanger.ac.uk/cosmic/mutation/overview'

snv_url = 'http://cancer.sanger.ac.uk/cosmic/gene/positive'

snv_params = {
    'src': 'gene',
    'coords': 'AA:AA',
    'fathmm': 'PATHOGENIC',
    'stat': ["1", "2"],
    'dr': '',
    'gd': '',
    'all_data': '',
    'start': 1,
    'export': 'json',
    'sEcho': 1,
    'iColumns': 19,
    'sColumns': '',
    'iDisplayStart': 0,
    'iDisplayLength': 800,
    'sn': 'upper_aerodigestive_tract',
    'bRegex': False,
    'iSortCol_0': 0,
    'sSortDir_0': 'asc',
    'iSortingCols': 1
}

reference_url = 'http://cancer.sanger.ac.uk/cosmic/references'

reference_params = {
    'q': 'MUTATION_STUD_REFERENCES',
    'repeat': 2,
    'export': 'json',
    'sEcho': 1,
    'iColumns': 7,
    'sColumns': '',
    'iDisplayStart': 0,
    'iDisplayLength': 1200,
    'sSearch': '',
    'bRegex': False,
    'iSortCol_0': 0,
    'sSortDir_0': 'asc',
    'iSortingCols': 1
}

headers = {
    'Accept-Encoding': 'gzip, deflate, sdch',
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_12_1) '
                  'AppleWebKit/537.36 (KHTML, like Gecko) Chrome/54.0.2840.98 Safari/537.36',
    'Referer': 'http://cancer.sanger.ac.uk/cosmic/gene/analysis',
    'Accept-Language': 'en,sv;q=0.8,en-US;q=0.6,zh-CN;q=0.4,zh;q=0.2',
    'Accept': 'application/json, text/javascript, */*; q=0.01',
    'X-Requested-With': 'XMLHttpRequest',
    'Connection': 'keep-alive'
    }

reference_filter = [
    'head and neck cancer',
    'head and neck cancers',
    'head and neck squamous cell carcinoma',
    'nasopharyngeal carcinoma',
    'oral cancer',
    'oral squamous cell carcinoma',
    'gingivo-buccal oral squamous cell carcinoma',
    'salivary gland tumours',
    'adenoid cystic carcinoma'
]

http = urllib3.PoolManager()


def decode(data):
    return data.decode('gbk', 'ignore').encode('utf-8')


def fetch_html(url, params=None, h=headers):
    if params is None:
        params = {}
    r = http.request('GET', url, fields=params, headers=h)
    return BeautifulSoup(decode(r.data), "html.parser")


def fetch_json(url, params=None, h=headers):
    if params is None:
        params = {}
    r = http.request('GET', url, fields=params, headers=h)
    return json.loads(r.data.decode('utf-8'))


def urlencode(params):
    return urllib.parse.urlencode(params, doseq=True)


def search_param(gen):
    params = search_params.copy()
    params['q'] = gen
    return params


def init_param(param, ln, gen_id, gen_seq):
    params = param.copy()
    if ln:
        params['ln'] = ln
    if gen_id:
        params['id'] = gen_id
    if gen_seq:
        params['seqlen'] = gen_seq
        params['end'] = gen_seq
    return params
